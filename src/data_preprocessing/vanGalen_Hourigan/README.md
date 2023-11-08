CORRECT THE NAME OF THE DIRECTORY. 
TITEL AND THE TEXT -> SAME STRUCTURE AS IN LASRY
# Studying Cell Type Communication Alterations in AML

To investigate alterations in cell-to-cell communication in AML, we constructed an integrated dataset containing the single-cell RNA sequencing (scRNAseq) profiles of bone marrow from healthy individuals ( [GSE120221 Oetjen et al., 2018](https://doi.org/10.1172/jci.insight.124928)) and AML patients at diagnosis ([GSE116256 vanGalen et al., 2019](https://doi.org/10.1016/j.cell.2019.01.031),). The raw read counts and annotation files were acquired from [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) and [GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221) datasets. Duplicated samples S1, Sk1, S2, Ck, and C2 were removed from the GSE120221 dataset. Additionally, the dataset GSE116256 had sample BM5-34p removed due to the absence of total bone marrow. Samples with less than 50% blasts were also excluded, resulting in 35 samples in the joint dataset.

## Data Preprocessing and Refinement

1. In the initial preprocessing step, we excluded genes that were missing in either of the datasets, resulting in 21,112 overlapping genes.
2. Genes with all-zero values in both datasets were removed, yielding 19,303 genes.
3. Cell IDs not present in cell annotation files were excluded, resulting in 87,333 cells.

## Quality Control and Filtering

1. To eliminate poor-quality cells, library size per cell was restricted to be between 1,000 and 30,000 reads.
2. The number of genes was required to be greater than 500.
3. Cells forming a cluster with high library size and low gene number were removed using a linear threshold of 3000*log10(library size + 1) - 10500.
4. This filtering process retained 78,119 cells.

## Cell Type Classification and Filtering

1. Small cell subtypes were grouped into 8 larger classes: HSPC, monocytes, DC, erythrocytes, megakaryocytes, B-cells, T-cells, and NK (Suppl Table N).
2. Mutation-bearing cells ("-like") and their healthy counterparts were merged into the same cell type categories.
3. Cell types were considered well-represented if they had at least 5 cells in each sample and were present in at least 30 samples.
4. NK and megakaryocytes did not meet these criteria and were excluded.
5. This filtering resulted in 6 cell types and 74,956 total cells.

## Gene Expression Filtering

1. Cell type pseudobulks were constructed by calculating the mean expression of a gene over cells within an individual cell type of a sample.
2. Cumulative pseudobulk counts were calculated for each gene by summing up the expression over the pseudobulks where it was detected.
3. Genes with a cumulative pseudobulk count greater than 1 were retained, resulting in 12,485 genes.

## Sample Filtering and Normalization

1. To ensure balanced cell type representation, samples with less than 5 cell types were removed.
2. This resulted in 33 samples and 74,583 total cells.
3. The 6 cell types displayed intrinsic differences in total RNA amounts.
4. Data normalization was performed using scran on each cell type separately, preserving natural differences.
5. Batch effect correction was carried out using scgen, and the batch-corrected counts were used for communication analysis.
