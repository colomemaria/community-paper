#!/bin/bash


# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    download_raw_data "Lasry"
    jupyter run /src/data_preprocessing/Lasry/1.preprocess_data.ipynb
    jupyter run /src/data_preprocessing/Lasry/2.filtering.ipynb
    jupyter run /src/data_preprocessing/Lasry/3.UMAP.ipynb
    
    jupyter run /src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/run_CPDB_Lasry.ipynb
    ./src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/runCPDB_Lasry.sh
elif [ "$1" == "--hourigan" ]; then
    download_raw_data "Hourigan"
else
    echo "Usage: $0 [--lasry | --hourigan]"
    exit 1
fi


