# Resource Usage: Custom Database

NicheNet requires 3 networks in order to generate a custom database and since we are building these networks via OmniPath, the first part of the script requires an internet access. Once the first part is done, we can run the Model Construction script and measure its resource usage with the below command. 

```bash
Rscript ModelConstruction.r outs/3_3/ & psrecord $(pgrep -f ModelConstruction.r) --interval 1 --log activity_CustomDB_NN.txt --plot plot__CustomDB_NN.png --include-children
```

The files that are required by CellPhoneDB can be found `/method_comparison/compare_algorithms/prepare_data/run_CellPhoneDB/build_customDB/CPDB_Custom` and measuring the resource usage can be done with the below command. 


```bash
cellphonedb database generate --user-interactions interactions.csv --user-interactions-only --user-protein prot_user.csv --user-gene gene_user.csv & psrecord $(pgrep -f cellphonedb) --interval 1 --log activity_CustomDB_CPDB.txt --plot plot__CustomDB_CPDB.png --include-children
```