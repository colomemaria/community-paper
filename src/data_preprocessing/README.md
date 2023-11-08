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
  - **Cell Filtering:** CELLS THAT ARE PRESENT INTHE COUNT MATRIX, BUT NOT IN THE CELL ANNOTATION FILE ARE FILTERED OUT
  - **Cell Type Filtering:** Excludes cell types with fewer than 5 cells per sample or present in fewer than 30 samples. THESE NUMBERS ARE DIFFERENT FOR DIFFERENT DATASETS, AREN'T THEY? I WOULD THUS NOT MENTION EXACT NUMBERS BUT JUST DESCRIBE WHAT IS DONE HERE.
  - **Gene Filtering:** Filters genes with low cumulative expression in constructed pseudo-bulks. AS FAR AS I REMEMBER, THE GENES ARE NOT FILTERED USING PSEUDOBULKS BUT ACTUAL CELLS -> PLEASE DOUBLE CHECK
  - **Sample Filtering:** Excludes samples with fewer than 5 cell types. AGAIN, DON'T USE THE EXCACT NUMBER SINCE IT MAY BE DIFFERENT FOR DIFFERENT DATASET, BUT EXPLAIN WHAT IS DONE IN GENERAL.
- **Output:** Processed data is saved in `.RData` format under `/results/data_preprocessing/$dataset/`.

### 3. Data Normalization
- **Notebook:** `./$dataset/3.normalization.ipynb`
- **Language:** R, using the Scran package LINK
- **Input:** Filtered data from the previous step.
- ANY OUTPUT YOU CREATE HERE? SAY THAT FOR SMILLIE, THIS IS THE OUTPUT THAT YOU USE FOR THE COMMUNITY. MENTION THAT LASRY DATA WAS ALREADY NORMALIZED, SO NO NORMALIZATION STEP WAS NEEDED. GIVE THE LINK TO THE ZENODO AGAIN, WHERE THE USE CAN DOWNLOAD THE FILES

### 4. Batch Correction
- **Notebook:** `./$dataset/4.1.batch_correction.ipynb`
- **Language:** Python, using the scgen library LINK
- **Input:** Normalized data from the previous step.
- **Warning:** This step may take a significant amount of time. It required nearly 14 hours on a system with 128GB RAM and 30 CPUs. For faster processing, consider using a GPU node.
- ANY OUTPUT YOU CREATE HERE? SAY THAT FOR CANGALEN-OETJEN, THIS IS THE OUTPUT THAT YOU USE FOR THE COMMUNITY. GIVE THE LINK TO THE ZENODO AGAIN, WHERE THE USE CAN DOWNLOAD THE FILES. PLEASE DOHBLE CHECK IF THE BATCH CORRECTION WAS ALSO DONE FOR SMILLIE AND LASRY, IF YES, SAY THAT IT WAS ONLY DONE FOR THE VISUALIZATION PURPOSES IN THE STEP 5, BUT NOT FOR THE COMMUNITY ANALYSIS.

### 5. Data Visualization
- **Notebook:** `./$dataset/4.2.visualization.ipynb`
- **Language:** R
- **Tasks:** This notebook visualizes the processed data, offering insights into the results of the preprocessing pipeline.


ADD A BRIEF SUMMARY OF WHAT FILES ARE PASSED TO COMMUNITY, WHAT IS THE STRUCTURE OF THESE FILES AND WHAT COLUMN NAMES ARE ESSENTIAL.
