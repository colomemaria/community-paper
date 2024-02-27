
# Integrated VanGalen-Oetjen Preprocessing - ScRNAseq Analysis of AML Bone Marrow
This repository contains the preprocessing pipeline and analysis steps for the integrated single-cell RNA sequencing (scRNAseq) profiles of bone marrow from healthy individuals ([GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221)) and AML (acute myeloid leukemia) patients at diagnosis ([GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256)). The raw data was sourced from [GSE120221 Oetjen et al., 2018](https://doi.org/10.1172/jci.insight.124928) and [GSE116256 vanGalen et al., 2019](https://doi.org/10.1016/j.cell.2019.01.031).

## Data Preprocessing Steps

1. **Data Source and Selection**:
   - The integrated dataset includes bone marrow samples from healthy individuals and AML patients.
   - Duplicated samples (`S1`, `Sk1`, `S2`, `Ck`, `C2`) and those not meeting specific criteria were removed. Additionally, the dataset GSE116256 with sample BM5-34p removed due to the absence of total bone marrow. Samples with less than 50% blasts were also excluded, resulting in 35 samples in the joint dataset.

2. **Initial Preprocessing and Gene Filtering**:
   - Excluded genes missing in both dataset, resulting in 21,112 overlapping genes.
   - Removed genes with all-zero values in both datasets, yielding 19,303 genes.

3. **Cell Filtering**:
   - Library size per cell was restricted between 1,000 and 30,000 reads.
   - Number of genes required to be greater than 500. This step retained 78,119 cells.

4. **Cell Type Classification and Filtering**:
   - Small cell subtypes grouped into 8 larger classes. Detailed information can be found in this [table.](https://github.com/colomemaria/community-paper/blob/main/data/cell_relabelling.csv)
   - Mutation-bearing cells and healthy counterparts were merged.
   - Cell types had to be present in at least 5 cells per sample and in at least 30 samples. This resulted in 6 cel types and 74.956 cells.

5. **Gene Filtering**:
   - Constructed cell type pseudobulks for gene filtering.
   - Genes with a cumulative pseudobulk count less than 1 were retained, resulting in 12,485 genes.

6. **Sample Filtering**:
   - Samples were assessed for balanced cell type representation and samples with less than 5 cell types were removed. The final dataset included 33 samples, a total of 74.583 cells.

7. **Data Normalization and Batch Correction:**
   - Data normalization performed using scran on each cell type separately.
   - Batch effect correction carried out using scgen.