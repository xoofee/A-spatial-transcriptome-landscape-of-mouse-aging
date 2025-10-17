# Spatial Transcriptomics Analysis Learning Roadmap

## Overview
This roadmap guides you through reimplementing the spatial transcriptomics analysis pipeline from "A spatial transcriptome landscape of mouse aging" (Cell 2024). The study examines 9 tissues (hippocampus, spinal cord, heart, lung, liver, small intestine, spleen, lymph node, testis) from male mice across 5 age points (2, 4, 13, 19, 25 months) using Stereo-seq technology. Key findings include: (1) increased organizational structure entropy (OSE) in aging tissues, (2) identification of senescence-sensitive spots (SSSs) as top 5% aging hotspots, (3) immunoglobulin G (IgG) accumulation as a hallmark of aging, (4) IgG-induced macrophage/microglia senescence, and (5) therapeutic potential of IgG reduction. The goal is to learn spatial analysis techniques by writing code from scratch, then referencing the original implementation when needed.

## Prerequisites
- Basic Python programming knowledge
- Understanding of single-cell RNA-seq concepts
- Familiarity with pandas, numpy, scipy basics
- Basic understanding of spatial transcriptomics (Stereo-seq)

## Learning Strategy
1. **Write code first** - Implement each step from scratch
2. **Reference original** - Check `./zhang-bin` code only when stuck
3. **Understand concepts** - Focus on the biological and statistical reasoning
4. **Test thoroughly** - Validate results against original paper

---

## Phase 1: Data Preprocessing & Quality Control (Week 1-2)

### 1.1 GEM File Processing
**Goal**: Convert raw GEM files to AnnData objects

**Tasks**:
- [ ] Write GEM file reader function using pandas
- [ ] Implement spatial binning algorithm (bin1 → bin50)
- [ ] Create expression matrix from spatial coordinates using scipy.sparse
- [ ] Generate AnnData spatial object with coordinates
- [ ] Add mitochondrial gene percentage calculation

**Key Learning Points**:
- GEM file format structure
- Spatial coordinate systems
- Binning algorithms for spatial data
- AnnData spatial object creation

**Reference**: `01_data_format/gem_to_seurat_obj.R`

### 1.2 Quality Control & Filtering
**Goal**: Implement comprehensive QC metrics

**Tasks**:
- [ ] Calculate QC metrics (n_genes, n_counts, pct_counts_mt) using scanpy
- [ ] Implement outlier detection (IQR method) with scipy.stats
- [ ] Create QC visualization plots using matplotlib/seaborn
- [ ] Apply filtering thresholds with pandas
- [ ] Generate QC reports using jupyter notebooks

**Key Learning Points**:
- Spatial transcriptomics QC metrics
- Outlier detection methods
- Visualization techniques for spatial data

---

## Phase 2: Spatial Analysis Fundamentals (Week 3-4)

### 2.1 Spatial Neighborhood Analysis
**Goal**: Understand spatial relationships between spots

**Tasks**:
- [ ] Implement neighbor detection algorithms using sklearn.neighbors
- [ ] Calculate spatial distance matrices with scipy.spatial
- [ ] Implement adjacency score analysis (K_z score) with permutation testing
- [ ] Create spatial graphs/networks using networkx
- [ ] Analyze spatial autocorrelation with scipy.stats
- [ ] Implement connected component analysis using scipy.ndimage

**Key Learning Points**:
- Spatial graph theory
- Neighborhood definitions
- Adjacency score calculation (K_z = (K_o - K_e) / K_sd)
- Permutation testing for spatial relationships
- Connected domain analysis

**Reference**: `02_spot_spatial_relationship/`

### 2.2 Spatial Segmentation with SLIC
**Goal**: Implement tissue segmentation using SLIC algorithm

**Tasks**:
- [ ] Implement SLIC (Simple Linear Iterative Clustering) algorithm from scratch
- [ ] Apply PCA for dimensionality reduction (PC1-PC15) with sklearn.decomposition
- [ ] Create spatial segments/superpixels with custom distance function
- [ ] Implement distance calculation: sqrt(m * d_spatial^2 + d_expression^2)
- [ ] Calculate segment properties and coverage with pandas
- [ ] Visualize segmentation results using matplotlib

**Key Learning Points**:
- SLIC algorithm implementation (iterative clustering)
- Combined spatial and expression distance metrics
- Superpixel generation with minimum spot requirements
- Coverage calculation using convex hull

**Reference**: `03_entropy/SLIC_segment_plot.R`

---

## Phase 3: Advanced Spatial Analysis (Week 5-6)

### 3.1 Organizational Structure Entropy (OSE) Analysis
**Goal**: Quantify spatial tissue organization disorder in aging

**Tasks**:
- [ ] Implement Organizational Structure Entropy (OSE) calculation
- [ ] Apply SLIC segmentation to create spatial segments
- [ ] Calculate Shannon entropy within each spatial segment: H = -Σ(p_ij * log2(p_ij))
- [ ] Implement Pielou's evenness: Pielou = entropy / log(condition)
- [ ] Normalize entropy by convex hull coverage: entropy / log(chull_range)
- [ ] Compare OSE scores between young (2M) and old (25M) mice
- [ ] Identify high OSE regions (R-OSEhigh) vs low OSE regions (R-OSElow)
- [ ] Create spatial OSE mapping visualizations

**Key Learning Points**:
- OSE as measure of tissue organization disorder
- SLIC-based spatial segmentation
- Shannon entropy calculation in spatial context
- Pielou's evenness index
- Convex hull coverage normalization
- Age-related spatial entropy changes
- Statistical comparison methods (Wilcoxon rank-sum test)

**Reference**: `03_entropy/tissue_SLIC_entropy.R`

### 3.2 Senescence-Sensitive Spots (SSSs) Analysis
**Goal**: Identify regions with heightened susceptibility to aging

**Tasks**:
- [ ] Implement aging score calculation using 5 well-established senescence gene sets
- [ ] Perform cell type-specific differential expression analysis (Wilcoxon test)
- [ ] Filter DEGs: p_val_adj < 0.05, |avg_log2FC| > 0.25, pct.1 > 0.1, pct.2 > 0.1
- [ ] Intersect aging DEGs with spot type-specific DEGs from each tissue
- [ ] Rank spots based on aging scores
- [ ] Identify top 5% spots as Senescence-Sensitive Spots (SSSs)
- [ ] Analyze distance-dependent aging score gradients (peaks at center, decreases outward)
- [ ] Implement SSS niche gradient analysis using Euclidean distance
- [ ] Analyze molecular gradients: humoral immune response, complement activation, SASP scores
- [ ] Correlate SSSs with OSE and cell identity scores

**Key Learning Points**:
- Senescence gene set integration (5 established sets)
- Cell type-specific DEG analysis
- Aging score calculation methodology
- Spatial hotspot identification (top 5%)
- Distance-dependent analysis
- Molecular gradient characterization
- SASP score correlation

**Reference**: `04_aging_related_DEGs/`

---

## Phase 4: Spatial Niche Analysis (Week 7-8)

### 4.1 Spatial Niche Analysis Around SSSs
**Goal**: Characterize microenvironments around senescence-sensitive spots

**Tasks**:
- [ ] Implement spatial niche detection around SSSs using Euclidean distance
- [ ] Calculate niche-specific markers using scanpy.tl.rank_genes_groups
- [ ] Implement niche gradient entropy analysis with resampling (R=100)
- [ ] Analyze niche composition and cell type distributions
- [ ] Compare niches across age groups and tissues
- [ ] Create niche visualization plots showing SSSs and surrounding areas
- [ ] Implement isoheight plots for SSS density visualization

**Key Learning Points**:
- Spatial niche concepts around hotspots
- Marker gene identification in niches
- Niche gradient entropy with resampling
- Compositional analysis of microenvironments
- Comparative spatial analysis across tissues
- Density-based visualization techniques

**Reference**: `05_SSS_nihce/`

### 4.2 Cell Identity Loss Analysis
**Goal**: Quantify cellular identity changes in aging tissues

**Tasks**:
- [ ] Implement cell identity scoring using AddModuleScore approach
- [ ] Calculate spatial cell type identity scores using top DEGs (top 20 genes)
- [ ] Filter identity genes: p_val_adj < 0.05, avg_log2FC > 0.25
- [ ] Analyze identity loss in R-OSEhigh vs R-OSElow regions
- [ ] Correlate identity scores with SSSs and aging
- [ ] Create identity visualization plots showing spatial patterns
- [ ] Implement temporal analysis across age groups (2M, 4M, 13M, 19M, 25M)

**Key Learning Points**:
- Cell identity quantification methods using module scoring
- Spatial context in cell typing
- Aging-related identity changes
- Correlation with OSE and SSSs
- Temporal analysis across multiple age points

**Reference**: `06_cell_identity/`

---

## Phase 5: Visualization & Integration (Week 9-10)

### 5.1 Advanced Visualization Techniques
**Goal**: Create comprehensive spatial visualizations

**Tasks**:
- [ ] Implement spatial overlap plots with multiple cell types
- [ ] Create isoheight plots for density visualization using stat_density_2d_filled
- [ ] Develop loess smoothing for temporal expression analysis
- [ ] Generate publication-quality multi-panel plots
- [ ] Implement spatial gradient visualizations
- [ ] Create interactive visualizations using plotly

**Key Learning Points**:
- Spatial overlap visualization techniques
- Density-based plotting methods
- Temporal smoothing algorithms
- Multi-panel plot design
- Gradient visualization
- Interactive spatial plots

**Reference**: `07_spatial_plot/`

### 5.2 Temporal Analysis and Multi-Tissue Integration
**Goal**: Integrate findings across 9 tissue types and temporal progression

**Tasks**:
- [ ] Implement temporal analysis across age groups (2M, 4M, 13M, 19M, 25M)
- [ ] Analyze IgG accumulation patterns across tissues and time
- [ ] Implement cross-tissue OSE comparison analysis
- [ ] Analyze tissue-specific aging patterns
- [ ] Compare SSSs across different tissues
- [ ] Create comprehensive multi-tissue visualizations
- [ ] Generate integrated aging landscape plots
- [ ] Implement loess smoothing for temporal trends

**Key Learning Points**:
- Multi-tissue integration methods
- Temporal analysis across 5 age points
- Cross-tissue comparison strategies
- Comprehensive aging landscape analysis
- Temporal smoothing techniques
- Statistical integration approaches

**Reference**: `08_bulk/`

---

## Key Methodologies from Original Code

### **SLIC Algorithm Implementation**
- **Distance Function**: `sqrt(m * d_spatial^2 + d_expression^2)` where m=3
- **Parameters**: window=20, min_spot=window^2*0.1, nPC=1:15
- **Iterative Process**: Update centers until convergence (res < 1)

### **Entropy Calculation**
- **Shannon Entropy**: `H = -Σ(p_ij * log2(p_ij))`
- **Pielou's Evenness**: `Pielou = entropy / log(condition)`
- **Coverage Normalization**: `entropy / log(chull_range)`

### **Adjacency Score Analysis**
- **Formula**: `K_z = (K_o - K_e) / K_sd`
- **Permutation Testing**: R=50 random samplings
- **Distance Threshold**: d=1 (nearest neighbors)

### **Cell Identity Scoring**
- **Method**: AddModuleScore with top 20 DEGs per cell type
- **Filtering**: p_val_adj < 0.05, avg_log2FC > 0.25
- **Temporal Analysis**: 2M, 4M, 13M, 19M, 25M age points

### **Key Findings from Paper**
- **OSE Increase**: Highest in hippocampus, spleen, lymph node, liver
- **SSS Cell Types**: Oligodendrocytes (hippocampus), neurons (spinal cord), cardiomyocytes (heart), plasmocytes (spleen/lymph node)
- **IgG Genes**: Igkc, Ighg1, Ighg2c, Ighg2b, Ighg3, Jchain most upregulated
- **Cross-Species**: IgG accumulation in human tissues (lymph node, liver, spleen, brain)
- **Therapeutic**: FcRn ASO reduces IgG and senescence markers
- **Mechanism**: IgG → Fcgr4 → macrophage/microglia senescence → tissue aging

---

## Phase 6: IgG Analysis & Therapeutic Interventions (Week 11-12)

### 6.1 Immunoglobulin G (IgG) Accumulation Analysis
**Goal**: Analyze IgG as a hallmark of aging

**Tasks**:
- [ ] Implement IgG expression analysis across tissues and age points
- [ ] Calculate IgG accumulation scores: Igkc, Ighg1, Ighg2c, Ighg2b, Ighg3, Jchain
- [ ] Analyze Ighigh spots (Igkc + heavy chain genes)
- [ ] Correlate IgG with SSSs and OSE scores
- [ ] Analyze IgG spatial patterns: co-localization, adjacency, distance from SSSs
- [ ] Validate IgG accumulation in plasma (ELISA)
- [ ] Cross-species validation: mouse (male/female) and human tissues
- [ ] Analyze IgG receptor expression: Fcgr4, Fcgrt

**Key Learning Points**:
- IgG as evolutionarily conserved aging biomarker
- Spatial IgG accumulation patterns
- Cross-species and cross-gender validation
- IgG receptor-mediated mechanisms

### 6.2 IgG-Induced Senescence Analysis
**Goal**: Investigate IgG as senescence-promoting factor

**Tasks**:
- [ ] Implement macrophage/microglia IgG treatment analysis
- [ ] Analyze senescence markers: SA-b-Gal, P21, MMTV, aggresomes
- [ ] Quantify inflammatory responses: NF-kB/P65, STAT1, iNOS, NO
- [ ] Implement Fcgr4 knockdown experiments
- [ ] Analyze human macrophage/microglia responses
- [ ] Validate heat-inactivated IgG controls
- [ ] Implement systemic IgG administration (6-month-old mice, 100 days)
- [ ] Analyze tissue IgG accumulation and senescence markers

**Key Learning Points**:
- IgG-induced cellular senescence mechanisms
- Macrophage/microglia senescence pathways
- Fcgr4 receptor-mediated signaling
- Systemic IgG effects on tissue aging

### 6.3 Therapeutic Intervention Analysis
**Goal**: Investigate IgG reduction as geroprotective strategy

**Tasks**:
- [ ] Implement FcRn (Fcgrt) ASO targeting
- [ ] Analyze IgG reduction in 19-month-old mice (40 days treatment)
- [ ] Quantify senescence marker reduction across tissues
- [ ] Compare with established geroprotective strategies (exercise, youthful fluids)
- [ ] Analyze IgG clearance mechanisms
- [ ] Evaluate therapeutic potential and limitations

**Key Learning Points**:
- FcRn-mediated IgG recycling
- ASO-based therapeutic interventions
- Geroprotective strategy validation
- Therapeutic target identification

---

## Phase 7: Advanced Topics & Validation (Week 13-14)

### 7.1 Method Validation
**Goal**: Validate your implementations against original results

**Tasks**:
- [ ] Compare results with original paper using pandas/scipy
- [ ] Implement statistical validation tests using scipy.stats
- [ ] Create reproducibility reports using jupyter notebooks
- [ ] Document differences and improvements using markdown
- [ ] Generate validation visualizations using matplotlib/seaborn

### 7.2 Advanced Spatial Analysis
**Goal**: Explore advanced spatial analysis techniques

**Tasks**:
- [ ] Implement spatial trajectory analysis using scanpy.tl.paga
- [ ] Apply spatial machine learning methods using sklearn
- [ ] Explore spatial network analysis using networkx
- [ ] Implement spatial clustering algorithms using sklearn.cluster
- [ ] Create advanced visualizations using matplotlib/seaborn/plotly

---

## Learning Resources

### Essential Python Packages
- `scanpy`: Single-cell and spatial analysis
- `pandas`, `numpy`: Data manipulation and numerical computing
- `scipy`: Scientific computing and statistics
- `matplotlib`, `seaborn`: Visualization
- `scikit-learn`: Machine learning algorithms
- `scikit-image`: Image processing
- `networkx`: Graph analysis
- `plotly`: Interactive visualizations
- `statsmodels`: Statistical modeling
- `jupyter`: Interactive notebooks

### Key Concepts to Master
1. **Spatial Data Structures**: Understanding coordinate systems
2. **Spatial Statistics**: Autocorrelation, clustering, segmentation
3. **Information Theory**: Entropy, mutual information
4. **Graph Theory**: Spatial networks, connectivity
5. **Machine Learning**: Clustering, classification in spatial context

### Validation Strategy
1. **Unit Testing**: Test each function independently
2. **Integration Testing**: Test complete pipelines
3. **Comparison Testing**: Compare with original results
4. **Visual Inspection**: Manually inspect key results
5. **Statistical Validation**: Use appropriate statistical tests

---

## Success Metrics

### Week 2: Basic Pipeline
- [ ] Successfully process GEM files to AnnData objects
- [ ] Generate QC plots and reports
- [ ] Understand spatial data structure

### Week 4: Spatial Analysis
- [ ] Implement neighborhood analysis
- [ ] Create spatial segmentation
- [ ] Understand spatial relationships

### Week 6: Advanced Analysis
- [ ] Calculate OSE scores and identify spatial organization changes
- [ ] Identify SSSs and analyze aging hotspots
- [ ] Understand senescence-associated spatial patterns

### Week 8: Niche Analysis
- [ ] Analyze spatial niches around SSSs
- [ ] Calculate cell identity loss patterns
- [ ] Understand microenvironment changes in aging

### Week 10: Visualization
- [ ] Create IgG spatial mapping visualizations
- [ ] Generate multi-tissue integration plots
- [ ] Understand comprehensive aging landscape visualization

### Week 12: IgG Analysis
- [ ] Analyze IgG accumulation patterns across tissues and ages
- [ ] Understand IgG-induced senescence mechanisms
- [ ] Implement therapeutic intervention analysis

### Week 14: Mastery
- [ ] Reproduce key paper results including IgG findings
- [ ] Understand all spatial analysis concepts and therapeutic implications
- [ ] Be able to apply methods to new datasets

---

## Troubleshooting Guide

### Common Issues
1. **Memory Problems**: Use sparse matrices, chunk processing
2. **Coordinate Systems**: Ensure consistent coordinate handling
3. **Statistical Testing**: Choose appropriate tests for spatial data
4. **Visualization**: Handle large datasets efficiently

### When to Reference Original Code
- Stuck on algorithm implementation
- Need to understand specific parameter choices
- Validation of intermediate results
- Understanding edge cases

### Getting Help
- Python documentation and tutorials
- Scanpy documentation and tutorials
- Spatial transcriptomics literature
- Online communities (Stack Overflow, Biostars)
- Original paper supplementary materials

---

## Final Notes

This roadmap is designed to build your understanding progressively. Each phase builds on the previous one, ensuring you understand both the technical implementation and the biological reasoning behind each analysis step.

Remember: The goal is not just to reproduce the code, but to understand the spatial analysis concepts deeply enough to apply them to new datasets and research questions.

Good luck with your spatial transcriptomics learning journey! 
