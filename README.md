# Community paper analysis

This repository contains a Python module and Jupyter Notebooks with R workflows for pre-processing of the single cell RNAseq datasets, which are used as showcases for cell-cell communication analysis with the `community` tool. The cell-cell communication analysis of these datasets can be found INSERT LINK TO THE COMMUNITY REPO HERE. YOU NEED TO MAKE CLEAR THAT THERE ARE TWO THINGS IN THIS REPO: 1. THE PRE-PROCESSIN OF THE DATASETS THAT WILL BE RUN WITH COMMUNITY -> LINK TO THE COMMUNITY REPO, 2. THE METHOD COMPARISON PART. TO BE HONSET, I DON'T SEE WHY WE NEED WHY WE NEED SEPARATELY THIS PART, AND RIGHT BELOW THE AIM / OVERVIEW PART. TO ME THE INFORMATION YOU HAVE IN BOTH THESE SECTIONS IS SIMILAR, SO I WOULD EVEN SUGGEST, YOU REMOVE THIS PARAGRAPH AND JUST SAY WAHT YOU WANT TO SAY IN THE AIM/OVERVIEW.

## Aim / Overveiw

YOU NEED TO MAKE CLEAR THAT THERE ARE TWO THINGS IN THIS REPO: 1. THE PRE-PROCESSIN OF THE DATASETS THAT WILL BE RUN WITH COMMUNITY -> LINK TO THE COMMUNITY REPO, 2. THE METHOD COMPARISON PART. 
The objective of this repository is to reproduce the analysis of the expression in cell-cell communication in disease state and identify potential differences in cell-cell communication in the diseased state. By doing so, the `community` tool will provide insights into the mechanisms of cell-cell communication and its potential role in the development of AML.


## Directory Structure

The repository is organized into two main directories: `/src` and `/results`. Below is the tree view of the `/src` directory which contains Jupyter Notebooks. The output from each section will be saved in the corresponding folder under the `/result/` directory.

![tree](https://imageupload.io/ib/JTjnBXgh9xdNnUa_1692782470.png)

IN THE FIGURE USE "DATA PRE-PROCESSING" INSTEAD OF "DATA CLEANSING". ALSO USE "COMPARE CELL COMMUNICATION RESULTS" INSTEAD OF "COMPARE RESULTS". YOU HAVE SEVERAL TYPES OF RESULTS IN THE METHOD COMPARISON, SO IF OYU JUST SAY "RESULTS", IT IS NOT PRECISE.


## Data pre-processing
WE DO NOT APPLY COMMUNITY IN THE RAW DATA SECTION. PLEASE CORRECT
We applied `community` tool on three published datasets.
    
1. Lasry et al, JOURNAL YEAR. The publication can be accessed via the following link: https://doi.org/10.1038/s43018-022-00480-0. The raw dataset can be downloaded by running `./download_raw_data.sh --lasry`. 
2. Similie et al, JOURNAL YEAR. This publication can be access via the following link: https://doi.org/10.1016/j.cell.2019.06.029. 
3. CORRECT THE NAME OF THE DATASET, MAKE SAME STYLE AS IN 1 AND 2: VanGalen_Hourigan, to study the alterations in cell-to-cell communication in AML, we constructed an integrated dataset containing the single-cell RNA sequencing (scRNAseq) profiles of bone marrow from healthy individuals ( [Oetjen et al., 2018](https://doi.org/10.1172/jci.insight.124928)) and AML patients at diagnosis ([vanGalen et al., 2019](https://doi.org/10.1016/j.cell.2019.01.031),). The raw read counts and annotation files were acquired from [GSE116256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256) and [GSE120221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221) datasets. The raw dataset can be downloaded by running `./download_raw_data.sh --hourigan`.

You can download both (YOU HAVE FOUR DATASETS -- WHICH DO YOU MEAN HERE??) datasets to corresponding directories by simply running

`make download-data`

**_NOTE:_** If you want to replicate the results and skip the preprocessing step of each dataset, you can download the pre-processed datasets from the following Zenodo links. 


### Supplementary tables -> IT BELONGS TO THE DATA PRE-PROCESSING SECTION, SO MAKE IT PART OF IT PLEASE

If need this section. YES, WE NEED THE CELL RELABELLING TABLE. PLEASE CONTACT FELIX AND UPLOAD THE TABLE HERE.


## Method comparison

DESCRIPTION OF THIS BRANCH HERE

## Resource usage

If this part is useful/neccessary, we need to check for this. YOU CAN SKIP THIS.

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
    
### Getting raw data -> TO THE DATA PRE-PROCESSING SECTION

IF YOU HAVE NOT YET GIVEN A LINK TO THE RAW DATA, DO IT IN THE VERY BEGINNING, PLEASE. ALSO IT WOULD BE GOOD IF YOU MADE IT CLEAR FOR EACH DOWNLOAD LINK YOU PUT WHETHER IT IS THE ORIGINAL RAW DATA FROM THE PAPERS OR THE OUTPUT FILES AFTER OUR PREPROCESSING. I AM A BIT CONFUSED HERE.
You can download the raw data to `/data/$dataset/raw_data` folder by running the below command. You can also visit the link here XXX_zenodo_link and download preprocessed files. 

- Download raw data into` /data/$dataset/raw_data` directory

    ```
    make download-data
    ```
