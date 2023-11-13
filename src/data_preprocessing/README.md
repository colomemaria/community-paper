## Preprocessing of Raw Data for Each Dataset

This section provides an overview of the preprocessing steps applied to each dataset.

### 1. Initial Data Cleaning and Annotation
- **Notebook:** `./$dataset/1.preprocess_data.ipynb`
- **Language:** R
- **Tasks:** This step handles initial cleaning and annotation of the raw data.

### 2. Data Filtering
- **Notebook:** `./$dataset/2.filtering.ipynb`
- **Language:** R
- **Tasks:** Further filtering of the processed data. Filtering includes:
  - **Cell Filtering:** Removes cells that appear in the count matrix but are absent in the cell annotation file.
  - **Cell Type Filtering:** Excludes cell types that do not meet the minimum prevalence criteria across samples, which may vary depending on the dataset characterics.
  - **Gene Filtering:** Filters out genes based on their cumulative expression accross cells, ensuring only genes with sufficient overall expression are retained for analyisis. 
  - **Sample Filtering:** Removes samples that fall below a certain threshold of cell type diversity, with specific riteria adjusted according to each dataset.
- **Output:** Processed data is saved in `.RData` format under `/results/data_preprocessing/$dataset/`. The Lasry dataset comes as normalized, therefore, the filtered count matrix from the this step is used. This pre-processed count matrix alongside with the ColData can be directly downloaded from [here](https://zenodo.org/records/7962808) 

### 3. Data Normalization
- **Notebook:** `./$dataset/3.normalization.ipynb`
- **Language:** R, using the [Scran package](https://bioconductor.org/packages/release/bioc/html/scran.html)
- **Input:** Filtered data from the previous step.
- **Output:** Normalized expression matrices. Smillie dataset does not need furhter pre-processing step and the output from this step is used as an input, as no further preprocessing is required. Downloadable files can be found on [Zenodo](https://zenodo.org/records/7962808)


### 4. Batch Correction
- **Notebook:** `./$dataset/4.1.batch_correction.ipynb`
- **Language:** Python, using the [scgen library](https://github.com/theislab/scgen)
- **Input:** Normalized data from the previous step.
- **Warning:** This step may take a significant amount of time. It required nearly 14 hours on a system with 128GB RAM and 30 CPUs. For faster processing, consider using a GPU node.
- **Output:** For the VanGalen-Oetjen dataset, this step is crucial, and the batch-corrected output is then used as input for the community tool. The pre-processed files are also stored on [Zenodo](https://zenodo.org/records/10013368)


### 5. Data Visualization
- **Notebook:** `./$dataset/4.2.visualization.ipynb`
- **Language:** R
- **Tasks:** This notebook visualizes the processed data, offering insights into the results of the preprocessing pipeline.

### Brief Summary of Structure of Input Files
**counts file:** normalized expression data frame containing gene symbols in the rows and cells in the columns.
**anno samples:** data frame of the sample annotation from all samples (rows are sample IDs, columns: must contain "sample_ID" and "case_or_control" columns).
**anno cells:** data frame of the cell annotation from all samples (rows are cell IDs, columns: must contain "cell_ID", "cell_type" and "sample_ID" columns).
