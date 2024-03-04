# libraries
library(nichenetr)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(igraph)
require(circlize)
library(R.utils)
library(data.table) #to read gz file
library(glue)


args <- commandArgs(trailingOnly = TRUE)


input_dir <- args[1]

#input_dir <- "../../../../../results/method_comparison/resoure_usage/Lasry/3_3/"

db_dir <- "../../../../compare_algorithms/prepare_data/run_NicheNet/build_customDB/NNET_Custom/"



# # load counts
# print("load counts")
# counts <- read.table(gzfile(paste0(path_in,"/counts_corr.csv.gz")
#                             )
#                      ,sep = ","
#                      ,row.names = 1
#                      ,header = TRUE
#                      )
# # load counts

counts <- fread(paste0(input_dir,"counts_norm.csv.gz"), header = TRUE,check.names=FALSE)
counts <- as.data.frame(counts)
rownames(counts) <- counts$V1
counts <- counts[,-1]
# head(str(counts))
# print(str(counts))

# load cell annotation
print("load cell annotation")
anno_cells <- read.table(paste0(input_dir,"anno_cells_norm.txt")
                         ,sep = "\t"
                         ,row.names = 1
                         ,header = TRUE
                         )
# print(str(anno_cells))

# load cell annotation
print("load cell annotation")
anno_samples <- read.table(paste0(input_dir,"anno_samples_norm.txt")
                         ,sep = "\t"
                         #,row.names = 2
                         ,header = TRUE
                         )
print(str(anno_cells))

#set rownames of annotation to cell_ids
rownames(anno_cells) <- anno_cells$cell_ID

#set colnames of counts to cell_ids
colnames(counts) <- rownames(anno_cells)

#anno_samples <- anno_samples[!anno_samples$region == "Non-inflamed",]
#anno_cells <- anno_cells[anno_cells$sample_ID %in% anno_samples$sample_ID,]
#counts <- counts[,colnames(counts) %in% rownames(anno_cells)]



srt=CreateSeuratObject(counts=counts, meta.data=anno_cells)

ligand_target_matrix = readRDS(paste0(db_dir,"ligand_target_matrixWithweights.rds"))
lr_network = readRDS(paste0(db_dir,"lig_rec_sources.rds"))
weighted_networks = readRDS(paste0(db_dir,"weighted_networksWithSourceWeights.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
# weighted_networks_lr = weighted_networks_lr %>% filter(from %in% lr_network$from & to %in% lr_network$to)






srt@meta.data$health_status %>% table()

srt@meta.data$cell_type %>% table()



Idents(srt) <- "cell_type"

cell_type_list <- unique(srt@meta.data$cell_type)
sample_list <- unique(srt@meta.data$sample_ID)

# # Create an empty list to store the differential expression tables for each cell type
DE_tables <- list()

# Iterate through each cell type in the cell_type_list
for (cell_type in cell_type_list){
    
    # Subset the Seurat object to the current cell type
    seurat_subset <- subset(srt, idents = cell_type)
    
    # Set the identity of the cells in the subset to the 'health_status' column
    seurat_subset <- SetIdent(seurat_subset, value = seurat_subset[["health_status"]])
    
    # Set the conditions of interest and reference for the differential expression analysis
    condition_oi <- "AML"
    condition_reference <- "healthy"
  
    # Call the function 'FindMarkers' to perform differential expression analysis and get the differential expression table
    DE_table <- FindMarkers(object = seurat_subset, ident.1 = condition_oi, 
                                    ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
    
    # Add the differential expression table for the current cell type to the 'DE_tables' list
    DE_tables[[cell_type]] <- DE_table
}


#DEGs <- readRDS("lasry_DEGs.rds")





# Create an empty list to store the expressed genes for each sample and cell type
expressed_genes_list <- list()

# Iterate through each sample in the sample_list
for (sample in sample_list){
    
    # Subset the Seurat object to the current sample
    sample_subset <- subset(x = srt, subset = sample_ID == sample)
    
    # Get a list of unique cell types in the current sample
    cell_types <- unique(sample_subset@meta.data$cell_type)
    
    # Iterate through each cell type in the current sample
    for (cell_type in cell_types){
        
        # Call the function 'get_expressed_genes' to get the genes that are expressed in the current cell type
        expressed_genes <- get_expressed_genes(cell_type, sample_subset, 0.10)
        
        # Convert the expressed genes to a vector and remove duplicates
        expressed_genes <- expressed_genes %>% unlist() %>% unique()
        
        # Add the list of expressed genes for the current sample and cell type to the 'expressed_genes_list'
        expressed_genes_list[[sample]][[cell_type]] <- expressed_genes
        
    }
}


sample_subset

length(expressed_genes_list[[sample]][[cell_type]])



# Define a function 'transform_interaction_df' that takes a data frame as input
transform_interaction_df <- function(df) {
  
  # Pivot the data frame to wide format, with 'interaction' as the ID column, 
  # 'sample' as the column names, and 'weight' as the values
  df_new <- df %>%
    pivot_wider(
      id_cols = interaction_ID,
      names_from = sample,
      values_from = c("weight")
    ) %>%
    
    # Reorder the columns to have 'interaction' as the first column
    select(interaction_ID, everything()) 
  
  # Return the pivoted data frame
  return(df_new)
}




# Get unique ligands and receptors from lr_network

ligands_in_db = lr_network %>% pull(from) %>% unique()
receptors_in_db = lr_network %>% pull(to) %>% unique()

# Initialize an empty list to store results
results=list()
ligand_activities_list <- list()
lr_network_df_large_full <- list()


for (receiver in cell_type_list){
    
    # Get DE genes for the receiver cell type
    #geneset_oi <- DE_tables[[receiver]] %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)


    if (nrow(DE_tables[[receiver]]) > 0) {
    
    # Get DE genes for the receiver cell type
    geneset_oi <- DE_tables[[receiver]] %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    } else {geneset_io <- "None"}


    
#     n_DE_genes <- length(geneset_oi)
    
    # Iterate over each sample
    for (sample in names(expressed_genes_list)){
        
        # Get expressed genes for the receiver cell type in the current sample
        expressed_genes_receiver <- expressed_genes_list[[sample]][[receiver]]
        
        
        # Get the intersection of DE genes and expressed genes for the receiver cell type in the current sample
        geneset_oi_sample_expressed <- geneset_oi[geneset_oi %in% expressed_genes_receiver]
        
        # Get the intersection of DE genes and ligand-target gene matrix rows
        geneset_oi_sample_expressed <- geneset_oi_sample_expressed[geneset_oi_sample_expressed %in% rownames(ligand_target_matrix)]
        
        # Get the intersection of expressed genes for the receiver cell type in the current sample and ligand-target gene matrix rows
        background_expressed_genes  <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)]
        
        # Iterate over each sender cell type
        for (sender in names(expressed_genes_list[[sample]])){
            
            
#             de_rec <- geneset_oi
#             de_lig <- DE_tables[[sender]] %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
            
            
            # Get expressed genes for the sender cell type in the current sample
            expressed_genes_sender <- expressed_genes_list[[sample]][[sender]]
            
            # Check if the sender cell type is not null for the receiver cell type
            if (!is.null(expressed_genes_list[[sample]][[receiver]])){
                
                # Get the intersection of expressed ligands and ligands in the lr_network
                expressed_ligands <- expressed_genes_sender[expressed_genes_sender %in% ligands_in_db]
                expressed_receptors <- expressed_genes_receiver[expressed_genes_receiver %in% receptors_in_db]
                
                
                
                # Get the potential ligands based on the expressed ligands and receptors
                potential_ligands = lr_network %>% 
                    filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
                        pull(from) %>% unique()
                
                
                
                potential_ligands = potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
                
                
                # Predict ligand activities with NicheNet
                ligand_activities = predict_ligand_activities(geneset = geneset_oi_sample_expressed, 
                                                              background_expressed_genes = background_expressed_genes, 
                                                              ligand_target_matrix = ligand_target_matrix, 
                                                              potential_ligands = potential_ligands)
                
                # Arrange ligand activities by pearson correlation coefficient and add a rank column
                ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
                
                
                # Get the potential ligands based on the ligand activities
                potential_ligands = ligand_activities %>% pull(test_ligand) %>% unique()
                
#                 potential_ligands = intersect(potential_ligands, de_lig)
                
                
#                 best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

                
                ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
                
                best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
                
                # Filter the lr_network to only include interactions between potential ligands and expressed receptors
                lr_network_expressed = lr_network %>% 
                    filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
                
                
                # Get the corresponding interactions from the weighted network
                lr_network_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & 
                                                                      to %in% expressed_receptors)
                
                # Add metadata to the network dataframe
                lr_network_df_large$sample <- sample
                lr_network_df_large$cell_type_l <- sender
                lr_network_df_large$cell_type_r <- receiver
                lr_network_df_large$pair <- paste(lr_network_df_large$from,lr_network_df_large$to, sep="_")
                lr_network_df_large$interaction_ID = paste0(lr_network_df_large$cell_type_l,":",
                                                         lr_network_df_large$from,"_",
                                                         lr_network_df_large$cell_type_r,":",
                                                         lr_network_df_large$to)
                
                
                df_new <- transform_interaction_df(lr_network_df_large)
                results[[sample]] <- rbind(results[[sample]],df_new)
                ligand_activities_list[[sample]][[sender]][[receiver]] <- ligand_activities
                
                lr_network_df_large_full[[sample]][[sender]][[receiver]] <- lr_network_df_large
#                 print(receiver)
#                 print(nrow(lr_network_df_large))
#                 write.csv(ligand_activities, paste("results6/",sample,"_",sender,"_",receiver,"_ligand_activities.csv", sep=""))
#                 write.csv(lr_network_df_large, paste("results6/",sample,"_",sender,"_",receiver,"_LR.csv", sep=""))
            }
        }
        
    }
    
}

