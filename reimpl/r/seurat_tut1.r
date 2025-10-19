
# install.packages(c("Seurat", "data.table", "Matrix", "ggplot2", "SeuratData"))

# devtools::install_github('satijalab/seurat-data')

#%%


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#%%

# InstallData("stxBrain")

# wget http://seurat.nygenome.org/src/contrib/stxBrain.SeuratData_0.1.2.tar.gz -OutFile stxBrain.SeuratData_0.1.2.tar.gz
# install.packages("./data/stxBrain.SeuratData_0.1.2.tar.gz", repos = NULL, type = "source")

brain <- LoadData("stxBrain", type = "anterior1")

# > head(brain@meta.data) class(brain@meta.data) = "data.frame"
#                    orig.ident nCount_Spatial nFeature_Spatial slice   region
# AAACAAGTATCTCCCA-1  anterior1          13069             4242     1 anterior
# AAACACCAATAACTGC-1  anterior1          37448             7860     1 anterior
# AAACAGAGCGACTCCT-1  anterior1          28475             6332     1 anterior
# AAACAGCTTTCAGAAG-1  anterior1          39718             7957     1 anterior
# AAACAGGGTCTATATT-1  anterior1          33392             7791     1 anterior
# AAACATGGTGAGAGGA-1  anterior1          20955             6291     1 anterior
# note: visium have barcode in the spot, it need not to generate x_1_y_1 like spot name

# > brain@assays           
# $Spatial
# Assay (v5) data with 31053 features for 2696 cells
# First 10 features:
#  Xkr4, Gm1992, Gm37381, Rp1, Sox17, Gm37323, Mrpl15, Lypla1, Gm37988,
# Tcea1
# Layers:
#  counts

# > rownames(brain@assays$Spatial)
#     [1] "Xkr4"             "Gm1992"           "Gm37381"
#     [4] "Rp1"              "Sox17"            "Gm37323"
#     [7] "Mrpl15"           "Lypla1"           "Gm37988"
#    [10] "Tcea1"            "Rgs20"            "Gm16041"
#    [13] "Atp6v1h"          "Oprk1"            "Npbwr1"
# ...

# > brain@images
# $anterior1
# Spatial coordinates for 2696 cells
# Default segmentation boundary: centroids
# Associated assay: Spatial
# Key: slice1_

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

#%%
coords <- GetTissueCoordinates(brain)

# Print first few rows
head(coords)

# > head(coords)
#                       x    y               cell
# AAACAAGTATCTCCCA-1 7475 8501 AAACAAGTATCTCCCA-1
# AAACACCAATAACTGC-1 8553 2788 AAACACCAATAACTGC-1
# AAACAGAGCGACTCCT-1 3164 7950 AAACAGAGCGACTCCT-1
# AAACAGCTTTCAGAAG-1 6637 2099 AAACAGCTTTCAGAAG-1
# AAACAGGGTCTATATT-1 7116 2375 AAACAGGGTCTATATT-1
# AAACATGGTGAGAGGA-1 8913 1480 AAACATGGTGAGAGGA-1


# make a scatter plot of the coordinates

ggplot(coords, aes(x = x, y = y)) +
  geom_point(size = 0.5, color = "blue") +
  scale_y_reverse() +       # flip Y to match tissue orientation
  coord_fixed() +           # keep X/Y scale equal
  theme_minimal() +
  labs(title = "Spatial Spot Coordinates", x = "X", y = "Y")
