# Lasry Preprocessing - ScRNAseq Analysis of Bone Marrow from Healthy and AML Individuals
This repository contains the preprocessing pipeline and analysis steps for the single-cell RNA sequencing (scRNAseq) data obtained from bone marrow samples of healthy and AML (acute myeloid leukemia) individuals. The raw data used in this study was sourced from publicly available dataset [Lasry, et al. 2022](https://www.nature.com/articles/s43018-022-00480-0) along with annotation files and downloaded from [GSE185381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381).

## Data Preprocessing Steps

1. **Data Source and Selection**:
   - In the AML cohort, we specifically included samples collected at diagnosis. 
   - Excluded duplicated samples from healthy individuals `healthy-4` and `healthy-5`. This resulted in a dataset of 14 samples, comprising both the healthy and AML cohorts, which serve as the input for the preprocessing pipeline. Additionally, genes with all-zero values were removed, resulting 31,843 genes.

2. **Cell Filtering:** Cells were filtered based on library size (between 1.100 and 30.00), and gene expression (more than 500 expressed genes). This filtering resulted in 57.420 cells remaining in the dataset.
3. **Cell Type Filtering**: Cell subtypes were grouped into 11 larger classes: HSPC, monocytes, granulocytes, DC, erythrocytes, megakaryocytes, perivascular cells, lymphoid progenitor cells (lymP), B-cells, T-cells, and NK cells. Detailed information can be found in this [table.](https://github.com/colomemaria/community-paper/blob/main/data/cell_relabelling.csv)

    Cell types with at least 5 cells per sample and presence in at least 12 samples were considered. Megakaryocytes, Perivascular cells, lymP cells, Natural killer, Dendritic cells, were excluded as they did not meet these criteria. 

    This resulted in 6 cell types and 56.662 total cells being retained.


4. **Gene Filtering**: Psuedo-bulk cells types were constructed (per sample) and filtered out genes that have too low cumulative expression in these pseudo-bulks. 15.770 genes left after this step.
5. **Sample Filtering:** Samples were assessed for balanced cell type representation. Those with fewer than 7 distinct cell types were excluded from the analysis. The final dataset included 13 samples, a total of 46.702 cells.

**Data Normalization:** The Lasry dataset was pre-normalized prior to our processing. Thus, no additional normalization steps were performed on this dataset as part of our preprocessing pipeline.