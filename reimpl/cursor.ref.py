

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%% Load GEM file
def load_gem_file(file_path) -> AnnData:
    """
    Load GEM file and convert to AnnData object
    
    Parameters:
    -----------
    file_path : str
        Path to GEM file (.gem or .gem.gz)
    
    Returns:
    --------
    adata : AnnData
        Spatial transcriptomics data
    """
    print(f"Loading GEM file: {file_path}")
    
    # Read GEM file
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
    else:
        df = pd.read_csv(file_path, sep='\t')
    
    print(f"GEM file shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    print(f"First few rows:")
    print(df.head())
    
    # Handle different column names
    if 'MIDCounts' in df.columns:
        df = df.rename(columns={'MIDCounts': 'MIDCount'})
    
    # Create expression matrix
    # Group by geneID and spatial coordinates, sum MIDCount
    expr_df = df.groupby(['geneID', 'x', 'y'])['MIDCount'].sum().reset_index()
    
    # Create sparse matrix
    # Get unique genes and spots
    genes = expr_df['geneID'].unique()
    spots = expr_df[['x', 'y']].drop_duplicates()
    spots['spot_id'] = range(len(spots))
    
    # Create mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}
    spot_to_idx = {(row['x'], row['y']): row['spot_id'] for _, row in spots.iterrows()}
    
    # Build sparse matrix
    row_indices = [gene_to_idx[row['geneID']] for _, row in expr_df.iterrows()]
    col_indices = [spot_to_idx[(row['x'], row['y'])] for _, row in expr_df.iterrows()]
    values = expr_df['MIDCount'].values
    
    # Create sparse matrix
    X = sp.csr_matrix((values, (row_indices, col_indices)), 
                      shape=(len(genes), len(spots)))
    
    # Create AnnData object
    adata = sc.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = [f"spot_{i}" for i in range(len(spots))]
    
    # Add spatial coordinates
    adata.obsm['spatial'] = spots[['x', 'y']].values
    
    # Add metadata
    adata.obs['x'] = spots['x'].values
    adata.obs['y'] = spots['y'].values
    
    print(f"AnnData shape: {adata.shape}")
    print(f"Spatial coordinates range: x=[{adata.obs['x'].min()}, {adata.obs['x'].max()}], y=[{adata.obs['y'].min()}, {adata.obs['y'].max()}]")
    
    return adata

#%% Quality Control Functions
def calculate_qc_metrics(adata):
    """
    Calculate quality control metrics
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data
    
    Returns:
    --------
    adata : AnnData
        Data with QC metrics added
    """
    print("Calculating QC metrics...")
    
    # Calculate basic metrics
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    # Calculate mitochondrial gene percentage
    mt_genes = adata.var_names.str.startswith('mt-')
    if mt_genes.sum() > 0:
        adata.obs['pct_counts_mt'] = (adata.X[:, mt_genes].sum(axis=1).A1 / 
                                     adata.obs['n_counts'].values * 100)
    else:
        adata.obs['pct_counts_mt'] = 0
    
    print(f"QC metrics calculated:")
    print(f"  n_genes: {adata.obs['n_genes'].describe()}")
    print(f"  n_counts: {adata.obs['n_counts'].describe()}")
    print(f"  pct_counts_mt: {adata.obs['pct_counts_mt'].describe()}")
    
    return adata

def detect_outliers_iqr(adata, metric='n_genes', factor=1.5):
    """
    Detect outliers using IQR method
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data
    metric : str
        Metric to use for outlier detection
    factor : float
        IQR factor for outlier detection
    
    Returns:
    --------
    adata : AnnData
        Data with outlier flags
    """
    print(f"Detecting outliers using IQR method for {metric}...")
    
    Q1 = adata.obs[metric].quantile(0.25)
    Q3 = adata.obs[metric].quantile(0.75)
    IQR = Q3 - Q1
    
    lower_bound = Q1 - factor * IQR
    upper_bound = Q3 + factor * IQR
    
    adata.obs[f'{metric}_outlier'] = (adata.obs[metric] < lower_bound) | (adata.obs[metric] > upper_bound)
    
    n_outliers = adata.obs[f'{metric}_outlier'].sum()
    print(f"  Found {n_outliers} outliers ({n_outliers/len(adata)*100:.1f}%)")
    print(f"  Bounds: [{lower_bound:.1f}, {upper_bound:.1f}]")
    
    return adata

def create_qc_plots(adata, save_path=None):
    """
    Create QC visualization plots
    
    Parameters:
    -----------
    adata : AnnData
        Spatial transcriptomics data
    save_path : str, optional
        Path to save plots
    """
    print("Creating QC plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Violin plots
    metrics = ['n_genes', 'n_counts', 'pct_counts_mt']
    for i, metric in enumerate(metrics):
        sns.violinplot(data=adata.obs, y=metric, ax=axes[0, i])
        axes[0, i].set_title(f'{metric} Distribution')
    
    # Spatial plots
    for i, metric in enumerate(metrics):
        scatter = axes[1, i].scatter(adata.obs['x'], adata.obs['y'], 
                                   c=adata.obs[metric], cmap='viridis', s=1)
        axes[1, i].set_title(f'{metric} Spatial Distribution')
        axes[1, i].set_xlabel('X coordinate')
        axes[1, i].set_ylabel('Y coordinate')
        plt.colorbar(scatter, ax=axes[1, i])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"QC plots saved to: {save_path}")
    
    plt.show()

#%% Main Analysis Pipeline
def main():
    """
    Main analysis pipeline
    """
    print("=== Spatial Transcriptomics Analysis Pipeline ===")
    print("Phase 1: Data Preprocessing & Quality Control")
    
    # File path
    gem_file = r"D:\work\bio\spatial\hallmark of aging\STT0000039\Analysis\STSA0000367\STTS0000720\Y-2M-Heart-1_1_filtered_bin50.gem.gz"
    
    # Check if file exists
    if not Path(gem_file).exists():
        print(f"Error: File not found: {gem_file}")
        return
    
    # Load GEM file
    adata = load_gem_file(gem_file)
    
    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Detect outliers
    adata = detect_outliers_iqr(adata, 'n_genes')
    adata = detect_outliers_iqr(adata, 'n_counts')
    adata = detect_outliers_iqr(adata, 'pct_counts_mt')
    
    # Create QC plots
    create_qc_plots(adata, 'qc_plots.png')
    
    # Save processed data
    adata.write('processed_spatial_data.h5ad')
    print("Processed data saved to: processed_spatial_data.h5ad")
    
    print("\n=== Phase 1 Complete ===")
    print("Next steps:")
    print("1. Apply quality filters")
    print("2. Implement spatial neighborhood analysis")
    print("3. Calculate organizational structure entropy (OSE)")
    print("4. Identify senescence-sensitive spots (SSSs)")
    
    return adata

#%% Run the analysis
if __name__ == "__main__":
    adata = main()
