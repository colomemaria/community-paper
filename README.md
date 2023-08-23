# Community paper analysis

This repository contains a Python module and Jupyter Notebooks with R workflows for processing and analyzing cell-cell communication data. The aim of this tool is to identify ligand-receptor interactions between cells and analyze potential differences in cell-cell communication in a diseased state, particularly in Acute Myeloid Leukemia (AML). The output of the analysis will include tables and figures.

## Aim / Overveiw

The objective of this repository is to reproduce the analysis of the expression in cell-cell communication in disease state and identify potential differences in cell-cell communication in the diseased state. By doing so, the `community` tool will provide insights into the mechanisms of cell-cell communication and its potential role in the development of AML.


## Raw data info
We applied `community` tool on three published datasets.
    
1. Lasry, the dataset associated with this research has undergone peer review and has been published in the journal Nature Cancer. The publication can be accessed via the following link: https://doi.org/10.1038/s43018-022-00480-0. The raw dataset can be downloaded by running `./download_raw_data.sh --lasry`. 
2. Similie, a published scRNA dataset of 48 biopsies taken from the colon of 12 healthy and 18 ulcerative colitis(UC) individuals. This publication can be access via the following link: https://doi.org/10.1016/j.cell.2019.06.029. 
3. VanGalen_Hourigan, to study the alterations in cell type to cell type communication in AML, we created an integrated dataset containing the scRNAseq of the bone marrow of healthy individuals (GSE116256 vanGalenREF, GSE120221 OetjenREF) and AML patients at diagnosis(REF). The raw dataset can be downloaded by running `./download_raw_data.sh --hourigan`. 

**_NOTE:_** If you want to replicate the results and skip the preprocessing step of each dataset, you can download the pre-processed datasets from the following Zenodo links. 


## Directory Structure

The repository is organized into two main directories: `/src` and `/results`. Below is the tree view of the `/src` directory which contains Jupyter Notebooks. The output from each section will be saved in the corresponding folder under the `/result/` directory.

![tree](https://imageupload.io/ib/JTjnBXgh9xdNnUa_1692782470.png)


## Supplementary tables

If need this section.


## Resource usage

If this part is useful/neccessary, we need to check for this.

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

    or
    
- Open:
    - [`src/1.preprocess_data.ipynb` Notebook](http://localhost:8888/notebooks/src/1.preprocess_data.ipynb) to run the demo workflow.
    
### Getting raw data

You can download the raw data to `/data/$dataset/raw_data` folder by running the below command. You can also visit the link here XXX_zenodo_link and download preprocessed files. 

- Download raw data into` /data/$dataset/raw_data` directory

    ```
    make download-data
    ```
