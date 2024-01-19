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
args <- commandArgs(trailingOnly = TRUE)


input_dir <- args[1]


sessionInfo()

data("LR_database")

# data("LR_database")
print(str(LR_database))



# # load counts
print("load counts")
counts <- fread(paste0(input_dir,"counts_corr.csv.gz"), header = TRUE)
counts <- as.data.frame(counts)
rownames(counts) <- counts$V1
counts <- counts[,-1]
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

# print(str(interactions))


# load("interactions.RData")



