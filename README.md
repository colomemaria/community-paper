# Community paper analysis

This repository contains a Python module and Jupyter Notebooks with R workflows for processing and analyzing cell-cell communication data. The aim of this tool is to identify ligand-receptor interactions between cells and analyze potential differences in cell-cell communication in a diseased state, particularly in Acute Myeloid Leukemia (AML). The output of the analysis will include tables and figures.

## Aim / Overveiw

The objective of this repository is to reproduce the analysis of the expression of one participant in cell-cell communication in AML and identify potential differences in cell-cell communication in the diseased state. By doing so, this tool will provide insights into the mechanisms of cell-cell communication and its potential role in the development of AML.


## Workflow

The raw data for this analysis is sourced from XXX, which is stored under the XXX_Zenodo_link. The first step in the analysis is performed using the R module in the Jupyter Notebook /src/1.preprocess_data.ipynb, which is responsible for cleaning and processing the raw data, as well as annotating the data using relevant files.

The processed data from the first step is then filtered using the Jupyter Notebook `/src/2.filtering.ipynb`. This R module filters out cells with low library size and low gene count, and the resulting processed data can be found in the /data/processed/ directory in .RData format. The figures generated from this step are stored in the `/results/filtered/figures/` directory.

The next step, performed using the Jupyter Notebook `/src/3.normalization.ipynb`, utilizes the Scran package in R to normalize the data. The input for this step is the filtered data from the previous step, which can be found in the `/data/filtered/` directory. The figures generated from this step are stored in the `/results/normalized/figures/` directory.

The Jupyter Notebook `src/4.1.batch_correction.ipynb` uses the scgen library in Python for batch correction. The input for this step is the normalized data from the previous step.

**Warning**
This notebook may take several hours to run, depending on the specifications of your computing system. It took almost 14 hours to run on 128GB RAM with 30 CPUs. 

`/src/4.2.visualization.ipynb` finally, last part of the workflow that is used to visualize the data before and after the batch correction. All the figures in here are save under `/results/visualization/figures/`

Finally, the Jupyter Notebook `/src/4.2.visualization.ipynb` is used to visualize the data before and after batch correction. All figures generated in this step are saved in the `/results/visualization/figures/` directory.

## Raw data info

Should we put some info about the raw data? 


## Supplementary tables

If need this section.


## Directory structure

The repository is organized into three main directories: `/src`, `/data`, and `/results`. The `/src` directory contains the Jupyter Notebooks that perform the analysis, `/data` directory contains the computed results from each analysis module, and `/results` directory is used to store the figures and tables generated by the respective modules.

The raw data should be downloaded into the /data/raw_data directory, which can be easily obtained by running the following command:

`make download-data`

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

You can download the raw data to `/data/raw_data` folder by running the below command. You can also visit the link here XXX_zenodo_link and download manually. 

- Download raw data into` /data/raw_data` directory

    ```
    make download-data
    ```
