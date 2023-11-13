# Aim / Overveiw

This repository is dedicated to the analysis of single-cell RNAseq datasets, facilitating a comprehensive understanding of cell-cell communication through two main components:

1. **Pre-Processing of Datasets:** Jupyter Notebooks are provided to prepare the datasets for subsequent analysis using the `community` tool. The pre-processed data is essential for accurate and effective cell-cell communication analysis.
2. **Method Comparison:** This section offers a comparative study of different tools (CellPhoneDB, NicheNet) in cell-cell communication analysis. It serves to highlight the advantages and use-cases of the `community` tool within the broader context of available analytical methods.

Additionally, this repository is designed in alignment with the FAIR (Findable, Accessible, Interoperable, and Reusable) principles to ensure the highest standards of open science. It serves as a resource for reproducing the results ~~presented in our associated paper~~, enabling researchers to validate and build upon our findings in the field of cell-cell communication

For detailed analysis using the `community` tool, refer to the [Community Repository](https://github.com/SoloveyMaria/community).

## Directory Structure

The repository is organized into two main directories: `/src` and `/results`. Below is the tree view of the `/src` directory which contains Jupyter Notebooks. The output from each section will be saved in the corresponding folder under the `/result/` directory.

![tree](https://imageupload.io/ib/2l26mryvfnbw3Ng_1699878309.png)


## Data pre-processing
We applied `community` tool on three published datasets.
    
1. Lasry et al, 2023. The publication can be accessed via the following link: https://doi.org/10.1038/s43018-022-00480-0. The raw dataset can be downloaded by running `./download_raw_data.sh --lasry`. 
2. Similie et al, 2019. This publication can be access via the following link: https://doi.org/10.1016/j.cell.2019.06.029. This published dataset is available for download from the controlled-access data repository, Broad DUOS, thus it requires a registration. https://portals.broadinstitute.org/single_cell/study/SCP259. 
3. Integrated van Galen et al. - Oetjen et al.: To study the alterations in cell-to-cell communication in AML, we constructed an integrated dataset containing the single-cell RNA sequencing (scRNAseq) profiles of bone marrow from healthy individuals ( [Oetjen et al., 2018](https://doi.org/10.1172/jci.insight.124928)) and AML patients at diagnosis ([vanGalen et al., 2019](https://doi.org/10.1016/j.cell.2019.01.031),). The  read counts and annotation files were acquired from [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) and [GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221) datasets. These both datasets can be downloaded by running `./download_raw_data.sh --hourigan`.


**_NOTE:_** If you want to replicate the results and skip the preprocessing step of each dataset, you can download the pre-processed datasets from the following Zenodo links. 

[Pre-processed data for Lasry et al, 2023](https://zenodo.org/records/7962808)

[Pre-processed data for Similie et al, 2019](https://zenodo.org/records/7962808)

[Pre-processed data for Integrated van Galen et al. - Oetjen et al](https://zenodo.org/records/10013368)

## Method comparison

This section is organized into four main parts, each designed to evaluate the `community` tool in relation to other established cell-cell communcation analysis tools. 

1. **Comparing Default Databases:** This part involves standardizing the databases to a common format to address the heterogeneity in gene space and ligand-receptor pair. Subsequently, we visualize and check the differences.
2. **Building a Unified Database and Running Algorithms:** We generate a customized database using `community` tool's database, ensuring each tool operates utilizing the identical database. 
3. **Comparing Results:** To facilitate a direct comparison, we first transform the output from each tool into a standardized format, overcoming the challenge of different results structures.
4. **Resource Usage:** The robustness and scability of the tools are put to the test through subsampling across a range from 6 to 32 scRNAseq samples, measuring the tools efficiency and adaptability to varying data sizes.


## How to run

- Clone the repo ```git clone https://github.com/colomemaria/community-paper.git``` and then cd into the directory ```cd community-paper/```

- If you don't have conda installed yet, install [conda](https://conda.io/miniconda.html) by running the command below

    ```
    make install-conda
    ```

- Create a conda environment named "community_paper" and install all necessary packages by using the following command:

    ```
    make create-env
    ```
- Launch [Jupyter](https://jupyter.org/) to access the notebooks to generate graphs

    ```
    make run-jupyter
    ```

- Go to [http://localhost:8888](http://localhost:8888) (a page should open automatically in your browser) 
    