# Lasry Preprocessing - ScRNAseq Analysis of Bone Marrow from Healthy and AML Individuals
This repository contains the preprocessing pipeline and analysis steps for the single-cell RNA sequencing (scRNAseq) data obtained from bone marrow samples of healthy and AML (acute myeloid leukemia) individuals. The raw data used in this study was sourced from publicly available dataset [Lasry, et al. 2022](https://www.nature.com/articles/s43018-022-00480-0) along with annotation files and downloaded from [GSE185381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381).

## Data Preprocessing Steps

1. **Data Source and Selection**:
   I DELETED THE FIRST TWO POINT WITH THE LINKS. YOU JUST GAVE THEM IN THE PARAGRAPH ABOVE -> NO NEED TO REPEAT IT HERE AGAIN.
   - In the AML cohort, we specifically included samples collected at diagnosis. 
   - Duplicated samples from healthy individuals #4 and #5 were excluded, resulting in a dataset of 14 samples. NAME THE EXACT IDS OF THE EXCLUDED SAMPLES, IF AVAILABLE. IS 14 SAMPLES THE HEALTHY COHORT OR THE HEALTHY + AML COHORT? IF IT IS THE FULL COHORT, PUT IT AS A SEPARATE SENTENCE BELOW, AND NOT TO THE POINT ABOUT THE HEALTHY

3. **Gene Filtering**:
   - In the initial pre-processing step, genes with all-zero values were removed.
   - This filtering process resulted in a dataset containing 31,843 genes.

4. **Cell Filtering**:
   - To ensure data quality, several filters were applied to cells, cell types, genes, and samples.
   - Cells with a library size between 1,100 and 30,000 reads were retained, and cells with greater than 500 expressed genes were kept.
   - After this step, 57,420 cells remained in the dataset.

5. **Cell Type Classification and Filtering**:
   - Cell subtypes were categorized into 11 larger classes, such as HSPC, monocytes, granulocytes, DC, erythrocytes, megakaryocytes, perivascular cells, lymP, B-cells, T-cells, and NK.
   - We merged mutation-bearing cells ("-like") and their healthy counterparts into the same cell type categories.
   - Cell types were considered well-represented if they had at least 5 cells in each sample and were captured in at least 12 samples.
   - Some cell types (megakaryocytes, perivascular cells, and lymP cells) did not meet these criteria and were excluded from the analysis.
   - As a result, 8 cell types and 56,662 total cells were retained in the dataset.

6. **Gene Expression Filtering**:
   - To filter out weakly expressed genes, cell type pseudobulks were constructed by calculating the mean expression of a gene over the cells in an individual cell type of a sample.
   - Cumulative pseudobulk counts were calculated for each gene by summing up the expression over the pseudobulks where it was detected.
   - Genes with a cumulative pseudobulk count greater than 0.25 were retained.
   - This filtering process resulted in 15,770 genes.

7. **Sample Filtering**:
   - To maintain a balanced cell type representation across all samples, a sample with less than 7 cell types was removed.
   - Ultimately, 13 samples and 46,702 total cells were included in the final dataset.

Please note that this repository is intended to document and share the preprocessing pipeline and data filtering steps undertaken in our study. It serves as a valuable resource for reproducibility and transparency in scRNAseq analyses.

