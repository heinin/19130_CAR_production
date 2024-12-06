#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 11/11/2024
# Description: ProjecTILs T cell annotation for the 19130 CAR products
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(ProjecTILs)
library(Seurat)
library(SignatuR)
library(UCell)

#==============================================================================#
# Helper functions
#==============================================================================#

source("/home/hnatri/Utilities/utilities.R")
source("/home/hnatri/13384_CART/CART_plot_functions.R")

#==============================================================================#
# Environment variables
#==============================================================================#

set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)

#==============================================================================#
# Import Seurat objects
#==============================================================================#

# Integrated data
seurat_data <- readRDS("/tgen_labs/banovich/pediatric_CAR-T/02_production/19130_CAR_production_integrated.rds")
seurat_data$cluster <- seurat_data$integratedSCTsnn_res.0.3

colnames(seurat_data@meta.data)

# Scaling and normalizing
DefaultAssay(seurat_data) <- "RNA"
#seurat_data <- NormalizeData(seurat_data)
#seurat_data <- ScaleData(seurat_data,
#                         vars.to.regress = c("percent.mt_RNA",
#                                             "percent.ribo_RNA",
#                                             "S.Score",
#                                             "G2M.Score"))
seurat_data <- FindVariableFeatures(seurat_data)

#saveRDS(seurat_data, "/tgen_labs/banovich/BCTCSF/TGFbR2_KO/TGFb_scRNAseq_integrated_seurat.rds")

query <- seurat_data
#DefaultAssay(query) <- "RNA"

#==============================================================================#
# Run ProjecTILs
#==============================================================================#

# Loading the reference atlas
ref <- load.reference.map()

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
p0 <- DimPlot(ref, label = T, cols = refCols) + coord_fixed(ratio=1)

# Expression of some marker genes across reference subtypes
markers <- c("Cd4", "Cd8a", "Ccr7", "Tcf7", "Pdcd1", "Havcr2", "Tox", "Izumo1r",
             "Cxcr6", "Xcl1", "Gzmb", "Gzmk", "Ifng", "Foxp3")
VlnPlot(ref, features = markers, stack = T, flip = T, assay = "RNA")

# Run Projection algorithm
query.projected <- Run.ProjecTILs(query, ref = ref)

# Visualize projection
p1 <- plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5) + coord_fixed(ratio=1)

# Plot the predicted composition of the query in terms of reference T cell subtypes
p2 <- plot.statepred.composition(ref, query.projected, metric = "Percent") + NoLegend()

# Compare gene expression
#plot.states.radar(ref, query = query.projected)

# Compare gene programs
programs <- GetSignature(SignatuR$Hs$Programs)
names(programs)

ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)
query.projected <- AddModuleScore_UCell(query.projected, features = programs, assay = "RNA",
                                        name = NULL)

#plot.states.radar(ref, query = query.projected, meta4radar = names(programs))

# The ProjecTILs.classifier function applies reference-projection but does not
# alter the current embeddings
querydata <- ProjecTILs.classifier(query = query, ref = ref)

palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(ref$functional.cluster), "NA")
p3 <- DimPlot(querydata, group.by = "functional.cluster", cols = palette) + coord_fixed(ratio=1)

(p1 + p2 + p3)

saveRDS(querydata, "/tgen_labs/banovich/pediatric_CAR-T/02_production/19130_CAR_production_integrated_ProjecTILs.rds")
#querydata <- readRDS("/tgen_labs/banovich/pediatric_CAR-T/02_production/19130_CAR_production_integrated_ProjecTILs.rds")

# Create barplot
palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd", "#777777")
names(palette) <- c(levels(ref$functional.cluster), "NA")
querydata$functional.cluster <- factor(querydata$functional.cluster, levels = levels(ref$functional.cluster))
# seurat_object, plot_var, group_var, group_levels, plot_levels, plot_colors, var_names, legend_title
barplot <- create_barplot(seurat_object = querydata,
                          plot_var = "functional.cluster",
                          group_var = "integratedSCTsnn_res.0.5",
                          group_levels = as.character(c(0, seq(1, max(unique(as.numeric(as.character(querydata$integratedSCTsnn_res.0.5))))))),
                          plot_levels = sort(unique(querydata$functional.cluster)),
                          plot_colors = palette,
                          var_names = c("% cells", "Cluster"),
                          legend_title = "T cell state")

barplot

