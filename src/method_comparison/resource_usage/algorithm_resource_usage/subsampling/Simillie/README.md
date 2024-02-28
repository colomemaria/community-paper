# Resource Usage: Algorithms, Lasry Simillie

To begin, we use the `split.py` script to subsample the data. This script takes the directory of the count matrix and splits it into smaller subsets under the given output directory.

```bash
./split.py ../../../../../data_preprocessing/Smillie/4.\ batch_correction/outs/ outs/
```

The following commands execute `community` on subsets of the Lasry dataset while simultaneously record resource usage with `psrecord` and we repeat the same procedues for each dataset and each tool. 


First subset using community

```bash
Rscript run_community.r outs/3_3/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Simillie_3_3_community.txt --plot plot_Smillie_3_3_community.png --include-children

```

Second subset using community

```bash
Rscript run_community.r outs/7_6/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Smillie_7_6_community.txt --plot plot_Smillie_7_6_community.png --include-children
```

Third subset using community

```bash
Rscript run_community.r outs/10_12/ & psrecord $(pgrep -f run_community.r) --interval 1 --log activity_Smillie_10_12_community.txt --plot plot_Smillie_7_6_community.png --include-children

```

First subset using CPDB and so on

```bash
Rscript PrepCPDB.r outs/3_3/ & psrecord $(pgrep -f PrepCPDB.r) --interval 1 --log activity_Smillie_3_3_cpdb.txt --plot plot_Smillie_3_3_cpdb.png --include-children```

Rscript PrepCPDB.r outs/7_6/ & psrecord $(pgrep -f PrepCPDB.r) --interval 1 --log activity_Smillie_7_6_cpdb.txt --plot plot_Smillie_3_3_cpdb.png --include-children```