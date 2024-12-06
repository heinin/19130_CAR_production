#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 08/28/2024
# Description: Inspecting CAR enrichment libraries
#==============================================================================#

#==============================================================================#
# Loading libraries and setting environment variables
#==============================================================================#

library(Seurat)
library(googlesheets4)

setwd("/home/hnatri/19130_CAR_production/")
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "umap"

source("/home/hnatri/SingleCellBestPractices/scripts/preprocessing_qc_module.R")

#==============================================================================#
# Import data
#==============================================================================#

# IND samples, CAR enrichment libraries only: F07628-CAR and F07629-CAR
path_list <- list("F07628" = "F07628-CAR_F07626-FB",
                  "F07629" = "F07629-CAR_F07627-FB")

# Without demultiplexing
# Importing data
seurat_list <- lapply(names(path_list), function(i){
  message(i)
  sample_10x_data <- Read10X(paste0("/scratch/emee/10x_fastq/Outs/New/", path_list[i], "/outs/filtered_feature_bc_matrix/"))
  sample_10x_data <- sample_10x_data$`Gene Expression`
  seurat_object <- CreateSeuratObject(counts = sample_10x_data)
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
})

names(seurat_list) <- names(path_list)

# Reads mapping to the construct
# Merge
seurat_list <- lapply(names(seurat_list), function(xx){
  seurat_obj <- seurat_list[[xx]]
  seurat_obj$FID <- xx
  
  seurat_obj
})
seurat_merged <- merge(x = seurat_list[[1]], y = seurat_list[1:length(seurat_list)])
seurat_merged <- JoinLayers(seurat_merged, assay = "RNA")

# The name of the construct is "IL13OP"
layer_data_RNA <- LayerData(seurat_merged,
                            assay = "RNA",
                            layer = "counts")
layer_data_RNA <- as.data.frame(t(as.matrix(layer_data_RNA)))

colnames(layer_data_RNA)[(length(colnames(layer_data_RNA))-10):length(layer_data_RNA)]
range(layer_data_RNA[,"IL13OP"])
range(layer_data_RNA)

layer_data_RNA$FID <- seurat_merged$FID

# What % of the library mapped to the construct?
layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)-1)))) %>%
  group_by(FID) %>%
  dplyr::summarise(CAR_counts = sum(IL13OP),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

# With demultiplexing
# Metadata
gs4_deauth()
TGFbRa2_KO_metadata  <- gs4_get("https://docs.google.com/spreadsheets/d/1T3oB6dS_rbrijBWKAlPUCOWNAFc2Zz1yHDsUJRSdpOg/edit?usp=sharing")
sheet_names(TGFbRa2_KO_metadata)
scRNAseq_samples <- read_sheet(TGFbRa2_KO_metadata, sheet = "scRNAseq")
scRNAseq_samples$Sample_Name_Run <- paste0(scRNAseq_samples$Sample_Name, "_",
                                           scRNAseq_samples$Run)

# Cell hashing antibody names
hash_antibodies <- unique(scRNAseq_samples$Hash_ID)

# Demultiplexing
seurat_data_list <- prep_seurat_list_multiplexed(metadata = scRNAseq_samples,
                                                 batch_ID = "Run",
                                                 cellRanger_path = "Enrichment_CellRanger_path",
                                                 cell_ID_prefix = "Sample_Name_Run",
                                                 CellHashing_Ab = "Hash_ID",
                                                 Hash_Abs = hash_antibodies,
                                                 pos_quant = 0.99)

seurat_data_list <- lapply(names(seurat_data_list), function(xx){
  obj <- seurat_data_list[[xx]]
  obj$FID <- xx
  
  obj
})

# Merging for plotting
seurat_data_merged <- merge(seurat_data_list[[1]], seurat_data_list[[2]])
seurat_data_merged <- JoinLayers(seurat_data_merged, assay = "RNA")

layer_data <- LayerData(seurat_data_merged,
                        assay = "RNA",
                        layer = "counts")
layer_data <- as.data.frame(t(layer_data))
layer_data$FID <- plyr::mapvalues(x = rownames(layer_data),
                                  from = rownames(seurat_data_merged@meta.data),
                                  to = as.character(seurat_data_merged@meta.data$FID))
layer_data$hash.ID <- plyr::mapvalues(x = rownames(layer_data),
                                      from = rownames(seurat_data_merged@meta.data),
                                      to = as.character(seurat_data_merged@meta.data$hash.ID))
layer_data$HTO_classification.global <- plyr::mapvalues(x = rownames(layer_data),
                                                        from = rownames(seurat_data_merged@meta.data),
                                                        to = as.character(seurat_data_merged@meta.data$HTO_classification.global))

table(seurat_data_merged$HTO_classification.global,
      seurat_data_merged$FID)

layer_data %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data)-3)))) %>%
  #dplyr::summarise(CAR_counts = sum(`CLTX-Construct`),
  #                 lib_size = sum(1:(ncol(layer_data)-2)),
  #                 pct_CAR = (CAR_counts/lib_size)*100) %>%
  group_by(FID) %>%
  dplyr::summarise(CAR_counts = sum(`IL13OP`),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

# What % of cells reach a threshold for CAR reads?
table(layer_data$FID)
nrow(layer_data)
layer_data %>% group_by(FID) %>% filter(`IL13OP` >= 1) %>% dplyr::summarise(nrow = n())
layer_data %>% group_by(FID) %>% filter(`IL13OP` >= 2) %>% dplyr::summarise(nrow = n())
layer_data %>% group_by(FID) %>% filter(`IL13OP` >= 3) %>% dplyr::summarise(nrow = n())

################################################################################

# Angela's data, no multiplexing
# Newly aligned, enrichment libraries only
# /scratch/emee/10x_fastq/Outs/New
# BCTCSF_0097_1_CS_Whole_C4_X5TTE_F04898_HVFC5DSX3
# BCTCSF_0097_1_CS_Whole_C4_X5TTE_F04902_HVFC5DSX3
path_list <- list(#"F04898" = "BCTCSF_0097_1_CS_Whole_C4_X5TTE_F04898_HVFC5DSX3",
                  "F04907" = "BCTCSF_0098_1_CS_Whole_C5_X5TTE_F04907_HVFC5DSX3",
                  #"F04902" = "BCTCSF_0097_1_CS_Whole_C4_X5TTE_F04902_HVFC5DSX3",
                  "F04900" = "BCTCSF_0099_1_CS_Whole_C2_X5TTE_F04900_HVFC5DSX3",
                  "F04906" = "BCTCSF_0097_1_CS_Whole_C6_X5TTE_F04906_HVFC5DSX3",
                  "F04904" = "BCTCSF_0099_1_CS_Whole_C2_X5TTE_F04904_HVFC5DSX3",
                  "F04884" = "BCTCSF_0097_1_CS_Whole_C7_X5TTE_F04884_HVFC5DSX3",
                  "F04908" = "BCTCSF_0099_1_CS_Whole_C5_X5TTE_F04908_HVFC5DSX3",
                  "F04899" = "BCTCSF_0098_1_CS_Whole_C3_X5TTE_F04899_HVFC5DSX3",
                  "F04901" = "BCTCSF_0103_1_PB_Whole_C1_X5TTE_F04901_HVFC5DSX3",
                  "F04903" = "BCTCSF_0098_1_CS_Whole_C3_X5TTE_F04903_HVFC5DSX3")
                  #"F04905" = "BCTCSF_0105_1_PB_Whole_C1_X5TTE_F04905_HVFC5DSX3")

sampletypes <- c("F04898" = "CSF",
                 "F04907" = "CSF",
                 "F04902" = "CSF",
                 "F04900" = "CSF",
                 "F04906" = "CSF",
                 "F04904" = "CSF",
                 "F04884" = "CSF",
                 "F04908" = "CSF",
                 "F04899" = "CSF",
                 "F04901" = "PBMC",
                 "F04903" = "CSF",
                 "F04905" = "Product")

# Importing data
seurat_list <- lapply(names(path_list), function(i){
  message(i)
  sample_10x_data <- Read10X(paste0("/scratch/emee/10x_fastq/Outs/New/", path_list[i], "/outs/filtered_feature_bc_matrix/"))
  seurat_object <- CreateSeuratObject(counts = sample_10x_data)
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
})

names(seurat_list) <- names(path_list)

# Reads mapping to the construct
# Merge
seurat_list <- lapply(names(seurat_list), function(xx){
  seurat_obj <- seurat_list[[xx]]
  seurat_obj$FID <- xx
  
  seurat_obj
  })
seurat_merged <- merge(x = seurat_list[[1]], y = seurat_list[1:length(seurat_list)])
seurat_merged <- JoinLayers(seurat_merged, assay = "RNA")

# The name of the construct is "IL13OP"
layer_data_RNA <- LayerData(seurat_merged,
                            assay = "RNA",
                            layer = "counts")
layer_data_RNA <- as.data.frame(t(as.matrix(layer_data_RNA)))

colnames(layer_data_RNA)[(length(colnames(layer_data_RNA))-10):length(layer_data_RNA)]
range(layer_data_RNA[,"IL13OP"])
range(layer_data_RNA)

layer_data_RNA$FID <- seurat_merged$FID

p1 <- layer_data_RNA %>% ggplot(aes(x = IL13OP)) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~FID)

p2 <- layer_data_RNA %>% ggplot(aes(x = log2(IL13OP))) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~FID)

p1 + p2

# What % of the library mapped to the construct?
layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)-1)))) %>%
  group_by(FID) %>%
  dplyr::summarise(CAR_counts = sum(IL13OP),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

layer_data_sum <- layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)))))

layer_data_sum[1:10, (ncol(layer_data_sum)-10):ncol(layer_data_sum)]

# What % of cells reach a threshold for CAR reads?
nrow(layer_data_RNA)
layer_data_RNA %>% filter(`IL13OP` >= 1) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 2) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 3) %>% nrow()

################################################################################

# Angela's samples, multiplexed
multipl_paths <- c("F05579" = "/scratch/emee/10x_fastq/Outs/New/BCTCSF_0114_1_PB_Whole_C1_X5TCR_F05579_HNJG2DSX5/",
                   "F04905" = "/scratch/emee/10x_fastq/Outs/New/BCTCSF_0105_1_PB_Whole_C1_X5TTE_F04905_HVFC5DSX3",
                   "F04901" = "/scratch/emee/10x_fastq/Outs/New/BCTCSF_0103_1_PB_Whole_C1_X5TTE_F04901_HVFC5DSX3")

# Importing data
seurat_list <- lapply(names(multipl_paths), function(i){
  message(i)
  sample_10x_data <- Read10X(paste0(multipl_paths[i], "/outs/filtered_feature_bc_matrix/"))
  seurat_object <- CreateSeuratObject(counts = sample_10x_data)
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt_RNA", assay = "RNA")
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]|^MRP[SL]", col.name = "percent.ribo_RNA", assay = "RNA")
})

names(seurat_list) <- names(multipl_paths)

# Reads mapping to the construct
# Merge
seurat_list <- lapply(names(seurat_list), function(xx){
  seurat_obj <- seurat_list[[xx]]
  seurat_obj$FID <- xx
  
  seurat_obj
})
seurat_merged <- merge(x = seurat_list[[1]], y = seurat_list[1:length(seurat_list)])
seurat_merged <- JoinLayers(seurat_merged, assay = "RNA")

# The name of the construct is "IL13OP"
layer_data_RNA <- LayerData(seurat_merged,
                            assay = "RNA",
                            layer = "counts")
layer_data_RNA <- as.data.frame(t(as.matrix(layer_data_RNA)))

colnames(layer_data_RNA)[(length(colnames(layer_data_RNA))-10):length(layer_data_RNA)]
range(layer_data_RNA[,"IL13OP"])
range(layer_data_RNA)

layer_data_RNA$FID <- seurat_merged$FID

p1 <- layer_data_RNA %>% ggplot(aes(x = IL13OP)) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~FID)

p2 <- layer_data_RNA %>% ggplot(aes(x = log2(IL13OP))) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(~FID)

p1 + p2

# What % of the library mapped to the construct?
layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)-1)))) %>%
  group_by(FID) %>%
  dplyr::summarise(CAR_counts = sum(IL13OP),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

layer_data_sum <- layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)))))

layer_data_sum[1:10, (ncol(layer_data_sum)-10):ncol(layer_data_sum)]

# What % of cells reach a threshold for CAR reads?
nrow(layer_data_RNA)
layer_data_RNA %>% filter(`IL13OP` >= 1) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 2) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 3) %>% nrow()

################################################################################

# Angela's data, GEX + CAR enrichment
product_pedi <- readRDS("/tgen_labs/banovich/pediatric_CAR-T/00_GEO_rds/product_T_cell_obj_2023.rds")

head(product_pedi@meta.data)
unique(product_pedi$FID_GEXFB)

# The name of the construct is "IL13OP"
layer_data_RNA <- LayerData(product_pedi,
                            assay = "RNA",
                            layer = "counts")
layer_data_RNA <- as.data.frame(t(layer_data_RNA))

range(layer_data_RNA[,"IL13OP"])
range(layer_data_RNA)

p1 <- layer_data_RNA %>% ggplot(aes(x = IL13OP)) +
  geom_histogram() +
  theme_classic()

p2 <- layer_data_RNA %>% ggplot(aes(x = log2(IL13OP))) +
  geom_histogram() +
  theme_classic()

p1 + p2

# What % of the library mapped to the construct?
layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA))))) %>%
  dplyr::summarise(CAR_counts = sum(IL13OP),
                   lib_size = sum(sum),
                   pct_CAR = (CAR_counts/lib_size)*100) %>%
  ungroup()

layer_data_sum <- layer_data_RNA %>%
  dplyr::mutate(sum = rowSums(across(1:(ncol(layer_data_RNA)))))

layer_data_sum[1:10, (ncol(layer_data_sum)-10):ncol(layer_data_sum)]

# What % of cells reach a threshold for CAR reads?
nrow(layer_data_RNA)
layer_data_RNA %>% filter(`IL13OP` >= 1) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 2) %>% nrow()
layer_data_RNA %>% filter(`IL13OP` >= 3) %>% nrow()

# CSF/PBMC
csf_pbmc_pedi <- readRDS("/tgen_labs/banovich/pediatric_CAR-T/00_GEO_rds/CSF_PBMC_all_cells_obj_2023.rds")

head(csf_pbmc_pedi@meta.data)

table(csf_pbmc_pedi$UPN_Cycle_Day_Sample_Type_Batch) %>%
  as.data.frame() %>%
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  RotatedAxis()

