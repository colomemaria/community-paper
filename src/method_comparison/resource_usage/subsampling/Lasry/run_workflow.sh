#!/bin/bash


# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    ./split.py ../../../../../results/data_preprocessing/Lasry/ ../../../../../results/method_comparison/resoure_usage/Lasry/
    
    
    Rscript lasry_community.r ../../../../../results/method_comparison/resoure_usage/Lasry/3_3/ & psrecord $(pgrep -f lasry_community.r) --interval 1 --log activity_Lasry_3_3_community.txt --plot plot_Lasry_3_3_community.png --include-children
    
    Rscript lasry_community.r ../../../../../results/method_comparison/resoure_usage/Lasry/7_6/ & psrecord $(pgrep -f lasry_community.r) --interval 1 --log activity_Lasry_7_6_community.txt --plot plot_Lasry_7_6_community.png --include-children
    
    
#    jupyter run /src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/run_CPDB_Lasry.ipynb
#    ./src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/runCPDB_Lasry.sh
elif [ "$1" == "--hourigan" ]; then
    download_raw_data "Hourigan"
else
    echo "Usage: $0 [--lasry | --hourigan]"
    exit 1
fi


