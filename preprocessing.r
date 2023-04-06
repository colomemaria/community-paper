library(dplyr)
library(Matrix)

################################################################################
#                              counts matrix                                   #
################################################################################
path_in <- paste0(getwd(),"/rdata")
print("input path is:")
print(path_in)

path_out <- paste0(getwd(), "/pdata")
print("output path is:")
print(path_out)
print("creating counts matrix")

file_names <- c("Epi", "Fib", "Imm")

gene_data <- lapply(file_names, function(file_name) {
  read.csv(paste0(path_in, "/", file_name, ".genes.tsv"), sep = "\t", header = FALSE)})
for (i in 1:length(gene_data)) {
  assign(paste0(file_names[i], "_genes"), gene_data[[i]])}

bar_data <- lapply(file_names, function(file_name) {
  t(read.csv(paste0(path_in, "/", file_name, ".barcodes2.tsv"), sep = "\t", header = FALSE))})
for (i in 1:length(bar_data)) {
  assign(paste0(file_names[i], "_bar"), bar_data[[i]])}

epi_all_genes <- union(union(epi_genes, fib_genes), imm_genes) # create rowname matrices
fib_all_genes <- union(union(fib_genes, imm_genes), epi_genes)
imm_all_genes <- union(union(imm_genes, epi_genes), fib_genes)

print("importing sorted matrices and creating counts submatrices")
all_genes <- lapply(file_names, function(file_name) {
  tmp_mat <- readMM(paste0(path_in, "/gene_sorted-", file_name, ".matrix.mtx")) %>%
    rbind(matrix(0, nrow = nrow(get(paste0(file_name,"_all_genes"))) - nrow(get(paste0(file_name,"_genes"))), 
                 ncol = ncol(get(paste0(file_name,"_bar"))))) %>%
    `colnames<-`(get(paste0(file_name,"_bar"))) %>%
    `rownames<-`(t(get(paste0(file_name,"_all_genes"))))
  assign(paste0(file_name, "_m"), tmp_mat, envir = .GlobalEnv)})

print("combining submatrices to counts matrix") #sorting submatrices and then combining them to counts matrix
counts <- cbind(epi_m[order(rownames(epi_m)),], fib_m[order(rownames(fib_m)),], imm_m[order(rownames(imm_m)),])

print("saving counts matrix")
save(counts, file = paste0(path_out,"/counts.RData"))
print("DONE")

################################################################################
#                              anno_samples                                    #
################################################################################

print("creating anno_samples")
print("importing metadata")
anno_samples <- read.csv(paste(file = paste0(path_in,"/all.meta2.txt")),sep = "\t",header = TRUE)
str(anno_samples)

print("correcting table")
anno_samples <- anno_samples[-1,] %>%
  rename(sample_ID = Sample,
         cell_ID = NAME,
         patient_ID = Subject,
         health_status = Health,
         cell_type_original = Cluster) %>%
  mutate(sample_ID = paste(sample_ID, health_status, sep = "_"))
str(anno_samples)

anno_samplex <- anno_samples[,-c(3,4,5,6,7)] #removing unnecessary rows
anno_samples <- anno_samples[,-c(2,3,4,7)]

anno_samples$case_or_control <- ifelse(grepl("healthy", anno_samples$health_status), "control", "case") #creating case/control

print("saving anno_samples as anno_samples.RData") #save anno_samples
save(anno_samples, file = "./pdata/anno_samples.RData")
print("DONE")

################################################################################
#                              anno_cells                                      #
################################################################################

print("creating anno_cells")
anno_cells <- read.csv(paste(file = './rdata/cell_relabelling.csv'),sep = ";")

print(str(anno_cells)) #rownames(cell_relabelling) <- cell_relabelling$cell_type_original
anno_cells <- anno_cells %>% right_join(anno_samplex, by=c("cell_type_original"), multiple = 'all')
anno_cells <- anno_cells[,-c(3,4,5)]
anno_cells <- anno_cells %>% relocate(cell_ID)

print("saving anno_cells as anno_cells.RData")
save(anno_cells, file = "./pdata/anno_cells.RData")
print("DONE")

################################################################################
#                              anno_genes                                      #
################################################################################

print("creating anno_genes")
load(file = "./rdata/LR_database.rda")
cat("str(LR_database)\n\n")

anno_genes <- data.frame(gene_symbol <- rownames(counts))

print("checking if genes are receptors or ligands")
anno_genes$inDB <- (anno_genes$gene_symbol %in% LR_database$Ligand)| (anno_genes$gene_symbol %in% LR_database$Receptor)
anno_genes$isLigand <- anno_genes$gene_symbol %in% LR_database$Ligand
anno_genes$isReceptor <- anno_genes$gene_symbol %in% LR_database$Receptor
str(anno_genes)

cat("total nr genes are", nrow(anno_genes),"\n")
cat("nr ligands in LR_database are", length(unique(LR_database$Ligand)),"\n")
cat("nr ligands in our data are",sum(anno_genes$isLigand),"\n")
cat("nr receptors in LR_database are",length(unique(LR_database$Receptor)),"\n")
cat("nr receptors in our data are",sum(anno_genes$isReceptor),'\n')

print("saving anno_genes as anno_genes.RData")
save(anno_genes, file = "./pdata/anno_genes.RData")
print("DONE")