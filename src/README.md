## Directory Structure

This section provides an overview of the directory structure for our GitHub repository, which is organized into two main sections: `data_processing` and `method_comparison`.

### 1. data_processing

The `data_processing` section focuses on processing the published raw data from Lasry and Similie. The data processing pipeline involves several steps such as:

- **Initial Cleaning:**  Cleaning and processing the raw data

- **Filtering:** LIST EXPLICITLY THE FOUR FILTERING STEPS. YOU ONLY MENTIONTWO OF THE HERE, AND THE MOST IMPORTANT ONES ARE NOT MENTIONED AT THE MOMENT. This step involves filtering out cell types that have too few cells, as well as identifying individual genes that may be specific to certain cell types under certain health conditions.

- **Normalization:** The data is normalized WITH SCRAN (LINK TO SCRAN PUBLICATION). The normalization is done within each cell type separately (-> PLEASE DOUBLE CHECK THAT IT IS TRUE FOR ALL DATASETS) in order to keep the cell-type intrinsic difference in the total RNA levels.

- **UMAP Visualization:** The processed data is visualized using the UMAP.
- YOU HAVE BATCH CORRECTION FOR VANGALEN-OETJEN. NAME THE METHOD YOU USED, ADD LINK TO THE PUBLICATION.

### 2. method_comparison

The `method_comparison` section is divided into three parts:

- **compare_databases:** This directory contains notebooks that compare the ORIGINNAL databases PROVIDED by each tool. 

- **compare_algorithms:** Here, in addition to running the tools, we construct a custom unified database for each tool, allowing us to compare the performance of different algorithms under standardized conditions. This ensures a fair and unbiased comparison.

- **compare_results:** The directory includes notebooks for visualizing and analyzing the output of each tool. This allows us to assess the results obtained from various methods and gain insights into their strengths and weaknesses.
