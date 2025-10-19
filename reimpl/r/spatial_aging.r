

# https://satijalab.org/seurat/articles/spatial_vignette

#%%

install.packages(c("Seurat", "data.table", "Matrix", "ggplot2"))


#%%
# ---------------------------------------------------------
# 🔧 Required packages
# ---------------------------------------------------------
library(Seurat)
library(data.table)
library(Matrix)
library(ggplot2)

# ---------------------------------------------------------
# 1️⃣ Load pre-binned GEM file (e.g., bin50.gem.gz)
# ---------------------------------------------------------
gem_file <- "D:/work/bio/spatial/hallmark of aging/STT0000039/Analysis/STSA0000400/STTS0000753/Y-2M-Liver-2_2_filtered_bin50.gem.gz"
cat("Loading bin50 GEM file...\n")

gem <- data.table::fread(gem_file)

# Expected columns:
# geneID, x, y, MIDCount, cell, ...
# Each (x, y) is already a 50x50µm bin coordinate
head(gem)

# ---------------------------------------------------------
# 2️⃣ Build expression matrix (genes × bins)
# ---------------------------------------------------------
expr_df <- gem[, .(count = sum(MIDCounts)), by = .(geneID, x, y)]

# Make spot_id from coordinates
expr_df[, spot_id := paste0("x", x, "_y", y)]

# Create sparse matrix
expr_mat <- sparseMatrix(
  i = as.integer(factor(expr_df$geneID)),
  j = as.integer(factor(expr_df$spot_id)),
  x = expr_df$count
)
rownames(expr_mat) <- levels(factor(expr_df$geneID))
colnames(expr_mat) <- levels(factor(expr_df$spot_id))

# ---------------------------------------------------------
# 3️⃣ Create Seurat object
# ---------------------------------------------------------
seu <- CreateSeuratObject(
  counts = expr_mat,
  project = "StereoSeq_bin50",
  assay = "Spatial"
)

# Add coordinates to metadata
coord_tbl <- unique(expr_df[, .(spot_id, x, y)])
coord_tbl <- coord_tbl[match(colnames(seu), coord_tbl$spot_id), ]
seu@meta.data$x <- coord_tbl$x
seu@meta.data$y <- coord_tbl$y

# ---------------------------------------------------------
# 4️⃣ QC filtering and normalization
# ---------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_Spatial > 100 & percent.mt < 20)

seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE)

# ---------------------------------------------------------
# 5️⃣ Dimensionality reduction, clustering, and visualization
# ---------------------------------------------------------
seu <- RunPCA(seu, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

# UMAP view
DimPlot(seu, reduction = "umap", label = TRUE) + ggtitle("UMAP Clusters")

# ---------------------------------------------------------
# 6️⃣ Spatial visualization
# ---------------------------------------------------------
meta <- seu@meta.data
p_spatial <- ggplot(meta, aes(x = x, y = y, color = factor(seu$seurat_clusters))) +
  geom_point(size = 0.8) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  ggtitle("Spatial Clusters (bin50)")

p_spatial

# ---------------------------------------------------------
# 7️⃣ Spatially variable features
# ---------------------------------------------------------
seu <- FindSpatiallyVariableFeatures(
  seu,
  assay = "SCT",
  features = VariableFeatures(seu)[1:2000],
  selection.method = "moransi"
)

top_features <- head(SpatiallyVariableFeatures(seu), 6)
FeaturePlot(seu, features = top_features)

# ---------------------------------------------------------
# 8️⃣ Save results
# ---------------------------------------------------------
saveRDS(seu, file = "StereoSeq_bin50_SeuratObj.rds")
