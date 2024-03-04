#!/bin/bash


# Check for the argument and call the function accordingly
if [ "$1" == "--lasry" ]; then
    ./split.py ../../../../../data_preprocessing/Lasry/2.filtering/outs/ outs/
    
    
    Rscript lasry_community.r outs/3_3/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_3_3_community.txt --plot plot_Lasry_3_3_community.png --include-children
    
    Rscript lasry_community.r outs/7_6/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_7_6_community.txt --plot plot_Lasry_7_6_community.png --include-children
    
    
#    jupyter run /src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/run_CPDB_Lasry.ipynb
#    ./src/method_comparison/compare_algorithms/run_CellPhoneDB/run_algorithm/runCPDB_Lasry.sh
else
    echo "Usage: $0 [--lasry | --hourigan]"
    exit 1
fi



