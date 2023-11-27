#!/bin/bash


# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    ./split.py 3 3 ../../../../results/method_comparison/resoure_usage/Lasry/3_3/ ../../../../results/data_preprocessing/Lasry/
    
    ./lasry_community.r
    
    
#    jupyter run /src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/run_CPDB_Lasry.ipynb
#    ./src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/runCPDB_Lasry.sh
elif [ "$1" == "--hourigan" ]; then
    download_raw_data "Hourigan"
else
    echo "Usage: $0 [--lasry | --hourigan]"
    exit 1
fi


