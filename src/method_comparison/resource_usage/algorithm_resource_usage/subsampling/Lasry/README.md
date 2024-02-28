# Resource Usage: Algorithms, Lasry Dataset

To begin, we use the `split.py` script to subsample the data. This script takes the directory of the count matrix and splits it into smaller subsets under the given output directory.

```bash
./split.py ../../../../../data_preprocessing/Lasry/2.filtering/outs/ outs/
```

The following commands execute `community` on subsets of the Lasry dataset while simultaneously record resource usage with `psrecord` and we repeat the same procedues for each dataset and each tool. 


First subset using community

```bash
Rscript run_community.r outs/3_3/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_3_3_community.txt --plot plot_Lasry_3_3_community.png --include-children

```

Second subset using community

```bash
Rscript run_community.r outs/7_6/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_7_6_community.txt --plot plot_Lasry_7_6_community.png --include-children
```

First subset using CPDB and so

```bash
Rscript PrepCPDB.r outs/3_3/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_3_3_cpdb.txt --plot plot_Lasry_3_3_cpdb.png --include-children```

Rscript PrepCPDB.r outs/7_6/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Lasry_7_6_cpdb.txt --plot plot_Lasry_7_6_cpdb.png --include-children```