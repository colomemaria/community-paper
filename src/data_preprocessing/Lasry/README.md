# Lasry Preprocessing - ScRNAseq Analysis of Bone Marrow from Healthy and AML Individuals
This repository contains the preprocessing pipeline and analysis steps for the single-cell RNA sequencing (scRNAseq) data obtained from bone marrow samples of healthy and AML (acute myeloid leukemia) individuals. The raw data used in this study was sourced from publicly available dataset [Lasry, et al. 2022](https://www.nature.com/articles/s43018-022-00480-0) along with annotation files and downloaded from [GSE185381](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185381).

## Data Preprocessing Steps

1. **Data Source and Selection**:
   I DELETED THE FIRST TWO POINT WITH THE LINKS. YOU JUST GAVE THEM IN THE PARAGRAPH ABOVE -> NO NEED TO REPEAT IT HERE AGAIN.
   - In the AML cohort, we specifically included samples collected at diagnosis. 
   - Duplicated samples from healthy individuals #4 and #5 were excluded, resulting in a dataset of 14 samples. NAME THE EXACT IDS OF THE EXCLUDED SAMPLES, IF AVAILABLE. IS 14 SAMPLES THE HEALTHY COHORT OR THE HEALTHY + AML COHORT? IF IT IS THE FULL COHORT, PUT IT AS A SEPARATE SENTENCE BELOW, AND NOT TO THE POINT ABOUT THE HEALTHY. SAY EXPLICITLY THAT IT IS THE NUMBER OF SAMPLES WICH IS THE INPUT TO THE PRE-PRCESSING PIPELINE, SUCH THAT PEOPLE DON'T THINK IT IS THE FINAL NUMBER

3. **Gene Filtering**:
   - In the initial pre-processing step, genes with all-zero values were removed.
   - This filtering process resulted in a dataset containing 31,843 genes. ALSO WEAKLY EXPRESSED GENES WERE REMOVED!! IS THIS NUMBER A FINAL NUMBER OR ONLY AFTER ZERO-VALUE GENE REMOVAL?? -> IF IT IS ONLY ABOUT REMOVING THE ZERO-VALUES GENES, I WOULD NOT MENTION IT HERE, AS IT IS CONFUSING. YOU HAVE THE GENE FILTERING LATER ON.

4. **Cell Filtering**:
   - To ensure data quality, several filters were applied to cells, cell types, genes, and samples. YOU ARE TALKING AXPLICITLY ABOUT THE CELL FILTERING HERE. WHY THIS SENTENSE ABOUT OTHER THREE TYPES OF FILTER? WHAT SENSE DOES IT MAKE HERE?
   - Cells with a library size between 1,100 and 30,000 reads were retained, and cells with greater than 500 expressed genes were kept.
   - After this step, 57,420 cells remained in the dataset. -> I LIKE THIS SENTENSE. IT IS VERY CLEAR. PLEASE USE THIS STYLE FOR THE OTHER FILTERS AS WELL. 

5. **Cell Type Filtering**:
   - Cell subtypes were categorized into 11 larger classes: HSPC, monocytes, granulocytes, DC, erythrocytes, megakaryocytes, perivascular cells, lymP, B-cells, T-cells, and NK. CELL RELABELLING CAN BE FOUND HERE (LINK TO THE CELL RELABELLING TABLE)
   - We merged mutation-bearing cells ("-like") and their healthy counterparts into the same cell type categories. ALL OTHER SENTENCES IN THIS SECTION ARE IN PASSIVE VIOCE -> DO THIS ONE IN THE PASSIVE VIOCE TOO. IN GENERAL, DON'T RANDOMLY SWITCH VIOCES IN THE TEXT. THINK FIRST WHAT SECTION SHOULD BE IN WHICH VOICE AND THEN STICK TO IT.
   - Cell types were considered well-represented if they had at least 5 cells in each sample and were captured in at least 12 samples. Megakaryocytes, perivascular cells, and lymP cells did not meet these criteria and were excluded from the analysis. 
   - As a result, 8 cell types and 56,662 total cells were retained in the dataset.

6. **Gene Filtering**: 
   - To filter out weakly expressed genes, cell type pseudobulks were constructed by calculating the mean expression of a gene over the cells in an individual cell type of a sample.
   - Cumulative pseudobulk counts were calculated for each gene by summing up the expression over the pseudobulks where it was detected.
   - Genes with a cumulative pseudobulk count greater than 0.25 were retained.
   - This filtering process resulted in 15,770 genes.

7. **Sample Filtering**:
   - To maintain a balanced cell type representation across all samples, a sample with less than 7 cell types was removed.
   - Ultimately, 13 samples and 46,702 total cells were included in the final dataset.

