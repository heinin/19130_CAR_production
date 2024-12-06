#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 10/25/2024
# Description: 19130 CAR production colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Colors and themes
#==============================================================================

# Colors for plotting
# Define colors for each level of categorical variables

# Cluster colors
product_clusters <- as.factor(c(0, seq(1:28)))

product_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(product_clusters))
names(product_cluster_col) <- levels(product_clusters)

# Cell type colors
gs4_deauth()
cluster_annot  <- gs4_get("https://docs.google.com/spreadsheets/d/1DstdK4vXjk7ZEeBalW_4CA6Vmd57BGENUx3iMDRLszk/edit?usp=sharing")
sheet_names(cluster_annot)
cluster_annot <- read_sheet(cluster_annot, sheet = "Cluster annotations")
head(cluster_annot)

product_celltype <- cluster_annot$cluster_name
product_celltype_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(product_celltype))
names(product_celltype_col) <- product_celltype

# Sample type
sample_type_col <- c("Mock" = "azure3",
                     "Pre_CAR" = "lightskyblue3",
                     "5050_CAR" = "darkslategray2",
                     "CAR" = "turquoise2")

day_col <- c("0" = "gray75",
             "1" = "darkolivegreen4",
             "4" = "darkolivegreen3",
             "7" = "darkolivegreen2",
             "14" = "darkolivegreen1")
