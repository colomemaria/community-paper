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
METHOD COMPARISON IS DONE ON THE XXX DATASET.

The `method_comparison` section is divided into three parts:
IT SHOULD HAVE THE FOLLOWING STRUCUTRE:
- **compare_databases:** This directory contains notebooks that compare the ORIGINNAL LIGAND-RECEPTOR databases PROVIDED by each tool. 
- **compare_cell_communication_results**
-- RUN XXX
-- RUN YYY
-- RUN ZZZ
-- COMPARE RESULTS

THE RUN XXX, RUN YYY, AND RUN ZZZ DIRECTORIES INCLUDE A STEP OF CONSTRUCTING A UNIFIED LIGAND-RECEPTOR DATABASE (BASED ON THE COMMUNITY DATABASE) TO COMPARE THE PERFORMANCE OF THE TOOLS UNDER STANDARDIZED CONDITIONS, AS WELL AS THE COMMUNICATION ANALYSIS BY EACH TOOL.

THE COMPARE RESULTS directory includes notebooks for visualizing and analyzing the output of each tool. 

