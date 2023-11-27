# # install community
# library(devtools)

# devtools::install_github("SoloveyMaria/community")

# libraries
library(community)
library(ggplot2)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(dendsort)
library(igraph)
require(circlize)
library(R.utils)
library(data.table) #to read gz file
library(Seurat)
library(nichenetr)

getwd()

input_dir <- "../../../../results/method_comparison/resoure_usage/Lasry/3_3/"

output_dir <- "output/"

sessionInfo()

data("LR_database")

# data("LR_database")
print(str(LR_database))



# # load counts
print("load counts")
counts <- readRDS("counts_clean/similie_counts.rda")
# counts <- as.data.frame(counts)
# rownames(counts) <- counts$gene_symbol
# counts <- counts[,-1]
print(str(counts))

# load cell annotation
print("load cell annotation")
anno_cells <- read.table(paste0(input_dir,"anno_cells_corr.txt")
                         ,sep = "\t"
                         ,row.names = 1
                         ,header = TRUE
                         )
print(str(anno_cells))

# load sample annotation
print("load sample annotation")
anno_samples <- read.table(paste0(input_dir,"anno_samples_corr.txt")
                           ,sep = "\t"
                           ,row.names = 1
                           ,header = TRUE
                           )
print(str(anno_samples))

colnames(counts) <- anno_cells$cell_ID
rownames(anno_cells) <- anno_cells$cell_ID





anno_cells <- anno_cells[anno_cells$sample_ID %in% anno_samples$sample_ID,]
counts <- counts[,colnames(counts) %in% rownames(anno_cells)]



seurat_obj=CreateSeuratObject(counts=counts, meta.data=anno_cells)

Idents(seurat_obj) <- "cell_type"

cell_type_list <- unique(seurat_obj@meta.data$cell_type)
sample_list <- unique(seurat_obj@meta.data$sample_ID)



# set threshold of the cell type size
threshold_celltype_size <- 6
print("threshold_celltype_size >")
print(threshold_celltype_size)

# set threshold of the minimum number of active cells
threshold_nr_active_cells <- 6
print("threshold_nr_active_cells >")
print(threshold_nr_active_cells)

# set threshold of expression
threshold_expr <- 0.1
print("threshold_expr >")
print(threshold_expr)

# Renaming the cell_ID.1 column in anno_cells to "cell_ID"
# colnames(anno_cells)[colnames(anno_cells) == "cell_ID.1"] <- "cell_ID"

# colnames(anno_cells)[colnames(anno_cells) == "cell_ID"] <- "cell_IDXX"



print("calculate communication")
interactions = calculate_communication(counts = counts
                                       ,anno_samples = anno_samples
                                       ,anno_cells = anno_cells
                                       ,threshold_celltype_size = threshold_celltype_size
                                       ,threshold_nr_active_cells = threshold_nr_active_cells
                                       ,threshold_expr = threshold_expr
                                       ,lrp_database = LR_database
                                       )

print("calculate general statistics")

interactions <- general_stat(comm_result = interactions
                                   ,verbose = FALSE#TRUE
)


threshold_log10_cum_weight <-  0.01
threshold_frac_samples_per_condition <-  0.6
threshold_log10_meanexpr_per_condition <- 0.02


print("filter weak interactions")

options(repr.plot.height = 10
       ,repr.plot.width = 16)
interactions <- filter_interactions(comm_result = interactions
                             ,threshold_frac_samples_per_condition = threshold_frac_samples_per_condition
                             ,threshold_log10_cum_weight = threshold_log10_cum_weight
                             ,threshold_log10_meanexpr_per_condition = threshold_log10_meanexpr_per_condition
)


threshold_log2FC <- 1
threshold_fdr <- 0.1


print("calculate differential communication")
interactions <- test_diff(comm_result = interactions
                          ,threshold_fdr = threshold_fdr
                          ,which_test = "t-test"
                          ,threshold_log2FC = threshold_log2FC
                          
                         )


# print(str(interactions))


# load("interactions.RData")



