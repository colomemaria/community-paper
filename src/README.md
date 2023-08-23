shit man how to start
## Directory Structure

This section provides an overview of the directory structure for our GitHub repository, which is organized into two main sections: `data_processing` and `method_comparison`.

### 1. data_processing

The `data_processing` section focuses on processing the published raw data from Lasry and Similie. The data processing pipeline involves several steps such as:

- **Initial Cleaning:**  Cleaning and processing the raw data

- **Filtering:** This step involves filtering out cell types that have too few cells, as well as identifying individual genes that may be specific to certain cell types under certain health conditions.

- **Normalization:** The data is normalized to ensure that different samples are comparable and that any biases are minimized.

- **UMAP Visualization:** The processed data is visualized using the UMAP, which helps in identifying patterns and relationships among the cells.

### 2. method_comparison

The `method_comparison` section is divided into three parts:

- **compare_databases:** This directory contains notebooks that compare the official databases used by each tool. By analyzing the databases, we can evaluate the quality and comprehensiveness of the data sources employed by different tools.

- **compare_algorithms:** Here, in addition to running the tools, we construct a custom unified database for each tool, allowing us to compare the performance of different algorithms under standardized conditions. This ensures a fair and unbiased comparison.

- **compare_results:** The directory includes notebooks for visualizing and analyzing the output of each tool. This allows us to assess the results obtained from various methods and gain insights into their strengths and weaknesses.
