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
  - **Cell Filtering:** Based on presence in the cell annotation file.
  - **Cell Type Filtering:** Excludes cell types with fewer than 5 cells per sample or present in fewer than 30 samples.
  - **Gene Filtering:** Filters genes with low cumulative expression in constructed pseudo-bulks.
  - **Sample Filtering:** Excludes samples with fewer than 5 cell types.
- **Output:** Processed data is saved in `.RData` format under `/results/data_preprocessing/$dataset/`.

### 3. Data Normalization
- **Notebook:** `./$dataset/3.normalization.ipynb`
- **Language:** R, using the Scran package
- **Input:** Filtered data from the previous step.

### 4. Batch Correction
- **Notebook:** `./$dataset/4.1.batch_correction.ipynb`
- **Language:** Python, using the scgen library
- **Input:** Normalized data from the previous step.
- **Warning:** This step may take a significant amount of time. It required nearly 14 hours on a system with 128GB RAM and 30 CPUs. For faster processing, consider using a GPU node.

### 5. Data Visualization
- **Notebook:** `./$dataset/4.2.visualization.ipynb`
- **Language:** R
- **Tasks:** This notebook visualizes the processed data, offering insights into the results of the preprocessing pipeline.
