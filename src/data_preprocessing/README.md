## Preprocessing of Raw Data for Each Dataset

In this section, we provide an overview of the preprocessing steps carried out for each dataset:

1. **Initial Data Cleaning and Annotation**:
   - The initial step is conducted using an R module within the Jupyter Notebook: `./$dataset/1.preprocess_data.ipynb`.
   - This module handles tasks such as cleaning and processing the raw data, as well as annotating the data using relevant files.

2. **Data Filtering**:
   - The processed data from the first step undergoes further filtering in the Jupyter Notebook: `./$dataset/2.filtering.ipynb`.
   - The filtering processes include:
        - Filter cells: by their presence in the cell annotation file
        - Filter cell types: Cell types which have too few cells (less than 5 cells) per sample or are present in too few samples (less than in 30 samples) is filtered out.
        - Filter genes: Individual genes might be specific to certain cell types in certain health conditions. We construct pseudo-bulk cell types (per sample) and filter out genes that have too low cumulative expression in these pseudo-bulks.
        - Filter samples: Samples with less than 5 cell types will be filtered out.
   - This R module is responsible for filtering out cells with low library size and low gene count.
   - The resulting processed data is stored in the `.RData` format within the `/results/data_preprocessing/$dataset/` directory.

3. **Data Normalization**:
   - The next stage, executed using the Jupyter Notebook: `./$dataset/3.normalization.ipynb`, employs the Scran package in R for data normalization.
   - The input for this stage is the filtered data obtained from the previous step, available in the `/results/data_preprocessing/$dataset/filtered/` directory.

4. **Batch Correction**:
   - The Jupyter Notebook: `./$dataset/4.1.batch_correction.ipynb` utilizes the Python-based scgen library for batch correction.
   - The input for this step is the normalized data from the previous normalization stage.

   **Warning:** Please note that running this notebook may take a considerable amount of time, depending on your computing system's specifications. It took nearly 14 hours on a system with 128GB RAM and 30 CPUs. For faster processing, consider running it on a GPU node.

5. **Data Visualization**:
   - Lastly, the Jupyter Notebook `4.2.visualization.ipynb` is employed to visualize the processed data, offering insights into the outcomes of the preprocessing pipeline.
