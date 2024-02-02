# notebook.output.textLineLimit=30

# libraries
library(nichenetr)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(igraph)
require(circlize)
library(R.utils)
library(data.table) #to read gz file

sessionInfo()

input_dir <- "../../../../../data_preprocessing/Lasry/2.filtering/outs/"
db_dir <- "../build_customDB/NNET_Custom/"
output_dir <- "outs/"



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
rownames(counts) <- counts$gene_symbol
counts <- counts[,-1]
head(str(counts))
print(str(counts))

anno_cells <- read.table(paste0(input_dir,"anno_cells_norm.txt")
                         ,sep = "\t"
                         ,header = TRUE
                         )

row.names(anno_cells) <- anno_cells$cell

colnames(counts) <- row.names(anno_cells)

head(colnames(counts))

head(rownames(anno_cells))

# get average expr 
avg_expression <- rowMeans(counts)
gene_names <- rownames(counts)
avg_expr_table <- data.frame(Gene = gene_names, OverallAverageExpression = avg_expression)

# get average expr within cohorts
case <- anno_cells$health_status == "AML"
case_idx <- rownames(anno_cells)[case]
avg_expression_case <- rowMeans(counts[,case_idx])


control <- anno_cells$health_status == "healthy"
control_idx <- rownames(anno_cells)[control]
avg_expression_control <- rowMeans(counts[,control_idx])

# average expr table
avg_expr_table$case_avg <- avg_expression_case
avg_expr_table$control_avg <- avg_expression_control
write_csv2(avg_expr_table, paste0(output_dir,"Lasry_avg_expr_table.csv"))

avg_expr_table



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


# total=0
# for (sample in sample_list){
    
#     for (cell_type in cell_type_list){
# #         print(sample)
# #         print(cell_type)
# #         print(length(expressed_genes_list[[sample]][[cell_type]]))
        
#         total=total+length(expressed_genes_list[[sample]][[cell_type]])
#     }
# }

# print(total)



sample_subset

length(expressed_genes_list[[sample]][[cell_type]])

# Create an empty list to store the differential expression tables for each cell type
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
    geneset_oi <- DE_tables[[receiver]] %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    
#     n_DE_genes <- length(geneset_oi)
    
    # Iterate over each sample
    for (sample in names(expressed_genes_list)){
        
        # Get expressed genes for the receiver cell type in the current sample
        expressed_genes_receiver <- expressed_genes_list[[sample]][[receiver]]
        
        
        # Get the intersection of DE genes and expressed genes for the receiver cell type in the current sample
        geneset_oi_sample_expressed <- geneset_oi[geneset_oi %in% expressed_genes_receiver]
        
        # Get the intersection of DE genes and ligand-target gene matrix rows
        geneset_oi_sample_expressed <- geneset_oi_sample_expressed[geneset_oi_sample_expressed %in% colnames(ligand_target_matrix)]
        
#         geneset_oi_sample_expressed <- geneset_oi_sample_expressed[geneset_oi_sample_expressed %in% lr_network$from]
        
        
        
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



matrix_result1 <- Reduce(
  
  # The `Reduce()` function takes two arguments: a function and a list.
  # In this case, the function is an anonymous function defined using the `function()` keyword.
  # This function takes two arguments `x` and `y` and performs a full join between them using the `full_join()` function from the `dplyr` package.
  # The `by = "interaction"` argument specifies that the join should be performed on the "interaction" column.
  function(x, y) full_join(x, y, by = "interaction_ID"), 
  
  # The second argument to the `Reduce()` function is a list called `results`.
  # This list contains data frames that need to be joined together.
  results
)

# Define a variable called `result` that will hold the output of the Reduce function
matrix_result <- Reduce(
  
  # The `Reduce()` function takes two arguments: a function and a list.
  # In this case, the function is an anonymous function defined using the `function()` keyword.
  # This function takes two arguments `x` and `y` and performs a full join between them using the `full_join()` function from the `dplyr` package.
  # The `by = "interaction"` argument specifies that the join should be performed on the "interaction" column.
  function(x, y) full_join(x, y, by = "interaction_ID"), 
  
  # The second argument to the `Reduce()` function is a list called `results`.
  # This list contains data frames that need to be joined together.
  results
)


# Example list of strings
strings <- matrix_result$interaction_ID

# Initialize empty vectors for each column
sender_celltype <- c()
sender_gene <- c()
receiver_celltype <- c()
receiver_gene <- c()

# Loop through each string and split it
for (string in strings) {
  parts <- strsplit(string, "_")
  
  # Split the sender part
  sender_parts <- strsplit(parts[[1]][1], ":")
  sender_celltype <- c(sender_celltype, sender_parts[[1]][1])
  sender_gene <- c(sender_gene, sender_parts[[1]][2])
  
  # Split the receiver part
  receiver_parts <- strsplit(parts[[1]][2], ":")
  receiver_celltype <- c(receiver_celltype, receiver_parts[[1]][1])
  receiver_gene <- c(receiver_gene, receiver_parts[[1]][2])
  interaction_ID <- c(string)
}

# Create a dataframe with the splitted values
df <- data.frame(
  sender_celltype = sender_celltype,
  sender_gene = sender_gene,
  receiver_celltype = receiver_celltype,
  receiver_gene = receiver_gene
)



# add interaction ID at the first column
df <- tibble::add_column(df, interaction_ID = matrix_result$interaction_ID, .before = "sender_celltype")



# Create interaction annotation df with log2FC values
for (row in 1:nrow(df)){
    each_row <- df[row,]
    sender_cell <- each_row$sender_celltype
    sender_gene <- each_row$sender_gene
    receiver_cell <- each_row$receiver_celltype
    receiver_gene <- each_row$receiver_gene
    
    DE_ligand_cell_type <- DE_tables[[sender_cell]]
    DE_receptor_cell_type <- DE_tables[[receiver_cell]]
    
    ligand_log2FC <- DE_ligand_cell_type[DE_ligand_cell_type$gene==sender_gene,]$avg_log2FC
    
    ligand_p_val_adj <- DE_ligand_cell_type[DE_ligand_cell_type$gene==sender_gene,]$p_val_adj
    
    receptor_log2FC <- DE_receptor_cell_type[DE_receptor_cell_type$gene==receiver_gene,]$avg_log2FC
    
    receptor_p_val_adj <- DE_receptor_cell_type[DE_receptor_cell_type$gene==receiver_gene,]$p_val_adj

    if (length(ligand_log2FC) == 0) {
        df[row,"ligand_log2FC"] <- NA
        df[row,"ligand_p_val_adj"] <- NA
    } else {
        df[row,"ligand_log2FC"] <- ligand_log2FC
        df[row,"ligand_p_val_adj"] <- ligand_p_val_adj
    }
    
    if (length(receptor_log2FC) == 0) {
        df[row,"receptor_log2FC"] <- NA
        df[row,"receptor_p_val_adj"] <- NA
    } else {
        df[row,"receptor_log2FC"] <- receptor_log2FC
        df[row,"receptor_p_val_adj"] <- receptor_p_val_adj
    }
}


#see where both components (sender/receiever) have log2FC value
df[complete.cases(df$ligand_log2FC, df$receptor_log2FC), ]

length(unique(lr_network$from))

for (gene in df$sender_gene){
    if (!gene %in% lr_network$from){
        print(gene)
    }
}



df[!complete.cases(df$ligand_log2FC, df$receptor_log2FC), ]



write.csv(matrix_result, paste0(output_dir,"NicheNet_scores.csv"))

write.csv(df, paste0(output_dir,"NicheNet_anno_interactions.csv"))

nn_interactions <- list()

nn_interactions$weights <- matrix_result
nn_interactions$anno_interactions <- df

# reassign directions for NN
nn_interactions$anno_interactions$direction_lig <- NA
nn_interactions$anno_interactions$direction_rec <- NA

nn_interactions$anno_interactions$direction_lig[!is.na(nn_interactions$anno_interactions$ligand_log2FC) & (nn_interactions$anno_interactions$ligand_log2FC < 0)] <- "down"
nn_interactions$anno_interactions$direction_lig[!is.na(nn_interactions$anno_interactions$ligand_log2FC) & (nn_interactions$anno_interactions$ligand_log2FC > 0)] <- "up"

nn_interactions$anno_interactions$direction_rec[!is.na(nn_interactions$anno_interactions$receptor_log2FC) & (nn_interactions$anno_interactions$receptor_log2FC < 0)] <- "down"
nn_interactions$anno_interactions$direction_rec[!is.na(nn_interactions$anno_interactions$receptor_log2FC) & (nn_interactions$anno_interactions$receptor_log2FC > 0)] <- "up"

nn_interactions$anno_interactions$direction_lig_rec <- paste(nn_interactions$anno_interactions$direction_lig
                                            ,nn_interactions$anno_interactions$direction_rec
                                            ,sep = "_")
idx_down <- (nn_interactions$anno_interactions$direction_lig_rec == "down_NA") | (
    nn_interactions$anno_interactions$direction_lig_rec == "NA_down") | (
    nn_interactions$anno_interactions$direction_lig_rec == "down_down")
idx_up <- (nn_interactions$anno_interactions$direction_lig_rec == "up_NA") | (
    nn_interactions$anno_interactions$direction_lig_rec == "NA_up") | (
    nn_interactions$anno_interactions$direction_lig_rec == "up_up")

nn_interactions$anno_interactions$direction <- NA
nn_interactions$anno_interactions$direction[idx_down] <- "down"
nn_interactions$anno_interactions$direction[idx_up] <- "up"
nn_interactions$anno_interactions$direction[!(idx_down | idx_up)] <- "ambigous"



save(nn_interactions, file="outs/nolog_nn_interactions.RData")


print("all good")
