#!/bin/bash

# get the the conda path
CONDA_PATH=$(dirname $(dirname $(which conda)))


# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    download_raw_data "Lasry"
    source $CONDA_PATH/bin/activate community_paper_latest
    jupyter run /src/data_preprocessing/Lasry/1.preprocess_data/1.preprocess_data.ipynb
    jupyter run /src/data_preprocessing/Lasry/2.filtering/2.filtering.ipynb
    jupyter run /src/data_preprocessing/Lasry/3.UMAP/3.UMAP.ipynb
    
    jupyter run /src/method_comparison/compare_algorithms/prepare_data/run_CellPhoneDB/build_customDB
    
    conda deactivate
    
    source $CONDA_PATH/bin/activate cellphonedb && cellphonedb database generate --user-interactions interactions.csv --user-interactions-only --user-protein prot_user.csv --user-gene gene_user.csv --result-path CPDB_Custom
    
    conda deactivate
    
    source $CONDA_PATH/bin/activate community_paper_latest
    
    jupyter run /src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/run_CPDB_Lasry.ipynb
    
    
    source $CONDA_PATH/bin/activate cellphonedb
elif [ "$1" == "--hourigan" ]; then
    download_raw_data "Hourigan"
else
    echo "Usage: $0 [--lasry | --hourigan]"
    exit 1
fi


