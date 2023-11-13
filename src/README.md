## Directory Structure

This section provides an overview of the directory structure for our GitHub repository, which is organized into two main sections: `data_processing` and `method_comparison`.

### 1. data_processing

The `data_processing` section focuses on processing the published data from Lasry et al., Similie et al., and Integrated van Galen et al. - Oetjen et al The data processing pipeline involves several steps such as: initial cleaning, filtering, normalization and batch correction (if applicable). For a detailed overview of our preprocessing steps, please go to (data_preprocessing directory)[./data_preprocessing].


### 2. method_comparison
METHOD COMPARISON IS DONE ON THE XXX DATASET. ->  SAY ON WHAT DATASET THE COMPARISON OF THE COMMUNICATION RESULTS WAS DONE AND ON WHAT DATASET THE RESURCE USAGE WAS DONE.

The `method_comparison` section is divided into three parts:
IT SHOULD HAVE THE FOLLOWING STRUCUTRE:
- **compare_databases:** Here, we explore and compare the original ligand-receptor databases as provided by each cell communication tool. The objective is to assess the variance between the databases and shared characteristics. 
- **compare_algorithms:** This directory is systemetically arranged into sub-dirs of each tool:
    * **run_run_NicheNet:** Contains the process of creating a standarized ligand-receptor database, derived from the `community` database, which is then used to analyze cell communication with NicheNet
    * **run_CellPhoneDB:** Foloows the same procedure as run_NicheNet, tailored for CellPhoneDB
    * **run_community**: running the community on the datasets. 

THE COMPARE RESULTS directory includes notebooks for visualizing and analyzing the output of each tool. 

