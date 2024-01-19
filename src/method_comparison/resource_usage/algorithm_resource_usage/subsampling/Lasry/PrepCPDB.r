# libraries
library(Seurat)
library(tidyverse)
library(igraph)
require(circlize)
library(R.utils)
library(data.table) #to read gz file
library(glue)


args <- commandArgs(trailingOnly = TRUE)


input_dir <- args[1]


input_dir <- input_dir
output_dir <- paste0(input_dir,"/CPDB/")
#final_out <- "../../../../../results/method_comparison/compare_algorithms/Lasry/CPDB/"
dir.create(file.path(output_dir))
print("here is the input_dir")
print(input_dir)

print("here is the output_dir")
print(output_dir)


# # load counts
# print("load counts")
# counts <- read.table(gzfile(paste0(path_in,"/counts_corr.csv.gz")
#                             )
#                      ,sep = ","
#                      ,row.names = 1
#                      ,header = TRUE
#                      )
# # load counts

counts <- fread(paste0(input_dir,"counts_corr.csv.gz"), header = TRUE,check.names=FALSE)
counts <- as.data.frame(counts)
rownames(counts) <- counts$gene_symbol
counts <- counts[,-1]
# head(str(counts))
# print(str(counts))

# load cell annotation
print("load cell annotation")
anno_cells <- read.table(paste0(input_dir,"anno_cells_corr.txt")
                         ,sep = "\t"
                         ,row.names = 1
                         ,header = TRUE
                         #,check.names=FALSE
                         )
# print(str(anno_cells))

#set rownames of annotation to cell_ids
rownames(anno_cells) <- anno_cells$cell_ID

#set colnames of counts to cell_ids
colnames(counts) <- rownames(anno_cells)

#create a Seurat object
srt=CreateSeuratObject(counts=counts, meta.data=anno_cells)


#peek into the number of cells for case/control
srt@meta.data$health_status %>% table()

#peek into the number of cell types
srt@meta.data$cell_type %>% table()

#set the indent to cell_type
Idents(srt) <- "cell_type"

# initialize empty vector for storing DEGs
DEGs <- readRDS("lasry_DEGs.rds")


meta <- anno_cells["cell_type"] %>% rownames_to_column("Cell")

# create a directory "samples_DEGs" to save the subsetted counts and annotation files. 
dir.create(file.path(output_dir, "samples_DEGs"))

# loop over each unique sample ID in the "sample_ID" column of the "anno_cells" data frame
for (sample in unique(anno_cells$sample_ID)) {
  
  # filter the annotation data frame to include only cells from the current sample
  anno_filtered <- filter(anno_cells, sample_ID == sample)
  
  # subset the expression counts matrix to the current sample
  subset_counts <- counts[, rownames(anno_filtered)]
  
  # subset the annotation data frame (required by CellPhoneDB)
  subset_meta <- anno_filtered["cell_type"] %>% rownames_to_column("Cell")
    
  # subset DEGs
  subset_DEGs <- DEGs %>% filter(cluster %in% unique(subset_meta$cell_type))
  
  # write the subsetted annotation data frame to a tab-separated value (TSV) file
  write.table(subset_meta, paste0(output_dir,"samples_DEGs/", sample, "_meta.tsv"), sep = '\t', quote = F, row.names = F)
  
  # write the subsetted counts matrix to a TSV file
  write.table(subset_counts, paste0(output_dir,"samples_DEGs/", sample, "_counts.tsv"), sep = '\t', quote = F)

  write.table(subset_DEGs, paste0(output_dir,"samples_DEGs/", sample, "_DEGs.tsv"), sep = '\t', quote = F)



}

output_dir <- paste0(output_dir, "/samples_DEGs/")

system(paste('conda run -n cpdb2 ./runCPDB.sh', output_dir))
