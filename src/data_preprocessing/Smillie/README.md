# Smillie Preprocessing - ScRNAseq Analysis of Colon Biopsies from Healthy and UC Individuals
This repository contains the preprocessing pipeline and analysis steps for the single-cell RNA sequencing (scRNAseq) data obtained from colon biopsies of healthy individuals and individuals with ulcerative colitis (UC). The raw data used in this study was sourced from the publicly available dataset [Smillie et al., 2019](https://doi.org/10.1016/j.cell.2019.06.029).

## Data Preprocessing Steps

1. **Data Source and Selection**:
   - The dataset consists of 48 biopsies collected from 12 healthy individuals and 18 UC patients. Biopsies from healthy patients included two separate areas of the colon, while UC patient samples encompassed both inflamed and non-inflamed tissue.
   - The original dataset contained 366,650 individual cells, with 365,492 cells utilized for preprocessing.

2. **Data Merging and Organization**:
   - Presorted epithelial, immune, and fibroblast submatrices were merged into a single counts matrix, with corresponding genes and cellular barcodes.

3. **Cell Filtering**:
   - Cells were filtered based on library size (between 1.100 and 30.00), and gene expression (more than 500 expressed genes). This step retained 160,707 cells out of 235,229.

4. **Cell Type Consolidation and Filtering**:
   - Original 50 cell types were consolidated into 19 broader categories. Detailed information can be found in this [table.](https://github.com/colomemaria/community-paper/blob/main/data/cell_relabelling.csv)
   - Cell types had to be present in at least 5 cells per sample and captured in 40 or more samples, resulting in 160,482 cells and 13 cell types.

5. **Gene Filtering**:
   - Cell type pseudobulks were constructed for filtering. Genes with a pseudobulk count greater than 3 were considered good quality, resulting in 13,861 genes out of 21,784.

6. **Sample Filtering**:
   - To address bias from duplicate samples, the one with the highest number of cell types was retained. For patient N51, sample B was retained due to irregular cell numbers in sample A.
   - This process yielded 28 samples (11 Healthy/17 Inflamed) and 100,000 cells for processing.

7. **Data Normalization**:
   - The filtered data underwent normalization by cell type using scran.