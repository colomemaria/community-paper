# libraries
library(Seurat)
library(tidyverse)
library(igraph)
require(circlize)
library(R.utils)
library(data.table) #to read gz file

input_dir <- "../../../../../results/data_preprocessing/Lasry/"
output_dir <- "../../../../../results/method_comparison/compare_algorithms/Lasry/CPDB/"
final_out <- "../../../../../results/method_comparison/compare_algorithms/Lasry/CPDB/"





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
                         ,check.names=FALSE
                         )
# print(str(anno_cells))

#set rownames of annotation to cell_ids
rownames(anno_cells) <- anno_cells$cell

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
DEGs <- c()

# iterate over each unique cell type 
for (cell in unique(srt@meta.data$cell_type)) {
  
  # subset Seurat object to only include cells of current cell type
  seurat_obj_receiver <- subset(srt, idents = cell)
  
  # set cell identity using the "health_status" feature
  seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["health_status"]])
  
  # specify the two conditions to compare
  condition_oi <- "AML"
  condition_reference <- "healthy" 
  
  # find differentially expressed genes between the two conditions
  DE_table_receiver <- FindMarkers(object = seurat_obj_receiver, 
                                   ident.1 = condition_oi, 
                                   ident.2 = condition_reference, 
                                   min.pct = 0.10) %>%
    # convert row names to a separate "gene" column
    rownames_to_column("gene")
  
  # add cell type information to the DEG table
  DE_table_receiver <- data.frame(cluster = cell, DE_table_receiver)
  
  # filter DEGs based on statistical significance and fold change threshold
  DE_table_receiver <- DE_table_receiver %>% 
    filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25)
  
  # print cell type and number of DEGs found
  print(cell)
  print(nrow(DE_table_receiver))
  
  # append DEGs to the vector of all DEGs
  DEGs <- rbind(DEGs, DE_table_receiver)
}


# write.table(DEGs, file =paste0(output_dir,"samples_DEGs/DEGs.tsv"), sep = '\t', quote = F, row.names = F)

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


#make sure cpdb is installed in the env.
run_CPDB <- './runCPDB.sh'

system('conda run -n cpdb ./runCPDB_Lasry.sh')

results_dir <- list.dirs(path = paste0(output_dir,"samples_DEGs/"), full.names = TRUE)

results_dir <- results_dir[grepl("_results", results_dir, fixed = TRUE)]

# Define a function called 'restructure_result' that takes one argument, 'cpdb_means'
restructure_result <- function(cpdb_means) {
  
  # Subset the columns of 'cpdb_means' that contain 'interacting_pair' or '|'
  cpdb_means <- cpdb_means[, grepl('interacting_pair|\\|', colnames(cpdb_means))]
  
  # Pivot the data to long format and split the 'interacting_pair' column into 'sending_protein' and 'receiving_protein' columns
  # Split the 'cell_types' column into 'sending_celltype' and 'receiving_celltype' columns
  # Unite the 'sending_celltype' and 'sending_protein' columns into a single column called 'sender'
  # Unite the 'receiving_celltype' and 'receiving_protein' columns into a single column called 'receiver'
  # Unite the 'sender' and 'receiver' columns into a single column called 'interacting_pairs'
  # Select the 'interacting_pairs' and 'value' columns
  conversion <- cpdb_means %>%
    pivot_longer(cols = -interacting_pair, names_to = "cell_types", values_to = "value") %>%
    separate(interacting_pair, c("sending_protein", "receiving_protein"), sep = "_") %>%
    separate(cell_types, c("sending_celltype", "receiving_celltype"), sep = "\\|") %>%
    unite(sender, c("sending_celltype", "sending_protein"), sep = ":", remove = FALSE) %>%
    unite(receiver, c("receiving_celltype", "receiving_protein"), sep = ":", remove = FALSE) %>%
    unite(interacting_pairs, c("sender", "receiver"), sep = "_", remove = FALSE) %>%
    select(interacting_pairs, value)
  
  # Return the processed data
  return(conversion)
}


results=list()
for (sample in results_dir){
    
    file <- paste0(sample,"/relevant_interactions.txt")
    
    sample_id <- basename(sample)
    sample_id <- strsplit(sample_id, '_')[[1]][1]
    
    
    if (file.exists(file)){
        
        cpdb_means <- read.csv(file, sep = "\t",  check.names = FALSE)
        
        
        sample_result <- restructure_result(cpdb_means)
        colnames(sample_result) <- c("interaction_ID",sample_id)
        results[[sample_id]] <- sample_result
        
    }
    
}

means=list()
for (sample in results_dir){
    
    file <- paste0(sample,"/means.txt")
    
    sample_id <- basename(sample)
    sample_id <- strsplit(sample_id, '_')[[1]][1]
    
    
    if (file.exists(file)){
        
        cpdb_means <- read.csv(file, sep = "\t",  check.names = FALSE)
        
        
        sample_result <- restructure_result(cpdb_means)
        colnames(sample_result) <- c("interaction_ID",sample_id)
        means[[sample_id]] <- sample_result
        
    }
    
}

# Define a variable called `result` that will hold the output of the Reduce function
means <- Reduce(
  
  # The `Reduce()` function takes two arguments: a function and a list.
  # In this case, the function is an anonymous function defined using the `function()` keyword.
  # This function takes two arguments `x` and `y` and performs a full join between them using the `full_join()` function from the `dplyr` package.
  # The `by = "interaction"` argument specifies that the join should be performed on the "interaction" column.
  function(x, y) full_join(x, y, by = "interaction_ID"), 
  
  # The second argument to the `Reduce()` function is a list called `results`.
  # This list contains data frames that need to be joined together.
  means
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

matrix_result[is.na(matrix_result)] <- 0

head(matrix_result)

str(matrix_result)



# str(matrix_result[rowSums(matrix_result[, -1] != 0, na.rm = TRUE) > 0, ])

str(matrix_result %>%
  filter(rowSums(. == 1) > 0))

matrix_result <- matrix_result[rowSums(matrix_result[, -1] != 0, na.rm = TRUE) > 0, ]

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
}

# Create a dataframe with the splitted values
df <- data.frame(
  sender_celltype = sender_celltype,
  sender_gene = sender_gene,
  receiver_celltype = receiver_celltype,
  receiver_gene = receiver_gene
)



# Create interaction annotation df with log2FC values
for (row in 1:nrow(df)){
    each_row <- df[row,]
    sender_cell <- each_row$sender_celltype
    sender_gene <- each_row$sender_gene
    receiver_cell <- each_row$receiver_celltype
    receiver_gene <- each_row$receiver_gene
    
    ligand_log2FC <- subset(DEGs, cluster == sender_cell & gene == sender_gene)$avg_log2FC
    
    ligand_p_val_adj <- subset(DEGs, cluster == sender_cell & gene == sender_gene)$p_val_adj
    
    receptor_log2FC <- subset(DEGs, cluster == receiver_cell & gene == receiver_gene)$avg_log2FC
    
    receptor_p_val_adj <- subset(DEGs, cluster == receiver_cell & gene == receiver_gene)$p_val_adj

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

df["interaction_ID"] <- paste0(df$sender_celltype, ":",df$sender_gene , "_", df$receiver_celltype, ":", df$receiver_gene)

df[is.na(df$ligand_log2FC) & is.na(df$receptor_log2FC),]

#subset means
means <- filter(means, interaction_ID %in% df$interaction_ID)

cols <- c("interaction_ID", grep("healthy", names(means), value = TRUE))
control_means <- means[, cols, drop = FALSE]

cols <- c("interaction_ID", grep("AML", names(means), value = TRUE))
case_means <- means[, cols, drop = FALSE]

numeric_cols <- control_means[, !names(control_means) %in% c("interaction_ID")]
control_means$row_means <- rowMeans(numeric_cols, na.rm = TRUE)

numeric_cols <- case_means[, !names(case_means) %in% c("interaction_ID")]
case_means$row_means <- rowMeans(numeric_cols, na.rm = TRUE)

df["mean_weigth_case"]=NA
df["mean_weight_control"]=NA



for (row in 1:nrow(df)){
    int_ID <- df[row,]$interaction_ID
    case_mean <- filter(case_means, interaction_ID==int_ID)$row_means
    control_mean <- filter(control_means, interaction_ID==int_ID)$row_means
    
    df[row,"mean_weigth_case"] <- case_mean
    df[row,"mean_weight_control"] <- control_mean
}

df["log2FC_weights"] = log2(df$mean_weigth_case/df$mean_weight_control)

names(df)

df <- df %>% select(interaction_ID, sender_celltype, sender_gene, receiver_celltype,receiver_gene,
             ligand_log2FC,ligand_p_val_adj,receptor_log2FC,receptor_p_val_adj,mean_weigth_case,
                   mean_weight_control,log2FC_weights)

threshold_log2FC <- 1

# Creating a new column named 'direction' based on the conditions

df$direction <- ifelse(df$log2FC_weights > threshold_log2FC, "up",
                             ifelse(df$log2FC_weights < threshold_log2FC, "down",
                                    ifelse(df$log2FC_weights == threshold_log2FC, "unchanged", NA)))


# The purpose of using a for loop in this code snippet is to handle the mismatched order of rows between 
# the means dataframe and the binary matrix dataframe.

# Create an empty list to store the multiplied rows
multiplied_rows <- list()

# Iterate over the rows
for (i in 1:nrow(means)) {
  interaction_ID <- means$interaction_ID[i]
  
  # Find the matching row in the 'significant' dataframe based on 'interaction_ID'
  matching_row <- matrix_result[matrix_result$interaction_ID == interaction_ID, ]
  
  # Perform element-wise multiplication
  multiplied_values <- means[i, -1] * matching_row[, -1]
  
  # Create a row with interaction_ID and multiplied values
  row <- c(interaction_ID, multiplied_values)
  
  # Add the row to the list
  multiplied_rows[[i]] <- row
}

# Convert the list of rows into a dataframe
multiplied_df <- do.call(rbind, multiplied_rows)

colnames(multiplied_df) <- c("interaction_ID", colnames(means)[-1])

multiplied_df <- as.data.frame(multiplied_df)

# Convert columns to double data type
multiplied_df <- as.data.frame(multiplied_df) %>%
  mutate(across(-interaction_ID, as.double))

multiplied_df$interaction_ID <- as.character(multiplied_df$interaction_ID)

multiplied_df[is.na(multiplied_df)] <- 0

library(community)

data(LR_database)

LR_DB <- LR_database



# LR_DB <- LR_DB %>% 
#         rename("Ligand" = "protein_name_a",
#                "Receptor" = "protein_name_b")

# LR_DB <- LR_DB[,-1]

df["pair"] <- paste0(df$sender_gene, "_", df$receiver_gene)
df["dup"] <- paste0(df$receiver_gene, "_", df$sender_gene)

# check if we have any duplicated swaps
df[df$pair %in% df$dup,]

# check if all the items that are not present in the original database exist as swapped pairs in 
#the original database.
identical(df[!df$pair %in% LR_DB$Pair.Name,]$interaction_ID,df[df$dup %in% LR_DB$Pair.Name,]$interaction_ID)

nrow(df)

fix_df = df[!df$pair %in% LR_DB$Pair.Name,]

df = df[df$pair %in% LR_DB$Pair.Name,]

nrow(df) + nrow(fix_df)

all(fix_df$dup %in% LR_DB$Pair.Name)

# create a pair column, makes it easier to check
multiplied_df$pair <- sapply(strsplit(multiplied_df$interaction_ID, "_"), function(x) {
  genes <- gsub(".*:", "", x)
  paste(genes, collapse = "_")
})

fix_multiplied_df = multiplied_df[!multiplied_df$pair %in% LR_DB$Pair.Name,]

multiplied_df = multiplied_df[multiplied_df$pair %in% LR_DB$Pair.Name,]

head(fix_multiplied_df)

# # Split values by underscore and swap
fix_multiplied_df$interaction_ID <- sapply(strsplit(fix_multiplied_df$interaction_ID, "_"), function(x) paste(rev(x), collapse = "_"))



colnames(fix_df)

new_df <- data.frame(
  interaction_ID = fix_df$interaction_ID,
  sender_celltype = fix_df$receiver_celltype,
  sender_gene = fix_df$receiver_gene,
  receiver_celltype = fix_df$sender_celltype,
  receiver_gene = fix_df$sender_gene,
  ligand_log2FC = fix_df$receptor_log2FC,
  ligand_p_val_adj = fix_df$receptor_log2FC,
  receptor_log2FC = fix_df$ligand_log2FC,
  receptor_p_val_adj = fix_df$ligand_p_val_adj,
  mean_weigth_case = fix_df$mean_weigth_case,
  mean_weight_control = fix_df$mean_weight_control,
  log2FC_weights = fix_df$log2FC_weights,
  direction = fix_df$direction
  )

new_df["interaction_ID"] <- paste0(new_df$sender_celltype, ":",new_df$sender_gene , "_", new_df$receiver_celltype, ":", new_df$receiver_gene)

multiplied_df <- rbind(fix_multiplied_df,multiplied_df)

df <- df[, !(names(df) %in% c("pair","dup"))]

df <- rbind(new_df,df)

all(multiplied_df$interaction_ID %in% df$interaction_ID)

multiplied_df <- multiplied_df[, !(names(multiplied_df) %in% c("pair","dup"))]

head(multiplied_df)



write.csv(multiplied_df, paste0(final_out,"CPDB_significant_weights.csv"))

write.csv(df, paste0(final_out,"CPDB_anno_interactions.csv"))







write.csv(matrix_result, paste0(final_out,"CPDB_weights.csv"))



anno <- read.csv(paste0(final_out,"CPDB_anno_interaction.csv"))

results <- read.csv(paste0(final_out,"CPDB_results.csv"))

threshold_log2FC <- 1

upregulated_anno <- anno[anno$log2>1,]

downregulated_anno <- anno[anno$log2<1,]

upregulated <- filter(results, interaction_ID %in% upregulated_anno$interaction_ID)

downregulated <- filter(results, interaction_ID %in% downregulated_anno$interaction_ID)



write.csv(upregulated_anno, paste0(final_out,"upregulated_anno.csv"))
write.csv(downregulated_anno, paste0(final_out,"downregulated_anno.csv"))
write.csv(upregulated, paste0(final_out,"upregulated.csv"))
write.csv(downregulated, paste0(final_out,"downregulated.csv"))



getwd()




