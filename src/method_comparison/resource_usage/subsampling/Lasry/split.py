#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scanpy as sc
import pandas as pd
import os
import shutil
import sys


n_case = int(sys.argv[1])
n_control = int(sys.argv[2])
path_out=sys.argv[3]
path_in = str(sys.argv[4])



data = sc.read_csv(path_in + "counts_corr.csv.gz")


anno_cells = pd.read_csv(path_in + "anno_cells_corr.txt", sep = "\t")


anno_samples = pd.read_csv(path_in + "anno_samples_corr.txt", sep = "\t")


case=anno_samples[anno_samples.case_or_control=="case"]



control=anno_samples[anno_samples.case_or_control=="control"]


data = data.T


anno_cells.index = anno_cells.cell_ID


data.obs = anno_cells.loc[data.obs.index]



case=case.sample(n_case)



control=control.sample(n_control)




ID_list = control.sample_ID.to_list() + case.sample_ID.to_list()



sub = data[data.obs.sample_ID.isin(ID_list)].copy()



sub.obs.sample_ID.value_counts()



sub = sub.T

path_out=sys.argv[3]

sub.to_df().to_csv(path_out + "counts_corr.csv.gz"
                   ,index=True
                   ,compression="gzip"
                   )




# export anno_cells_corr
print("save anno_cells_corr.txt")
sub.var.to_csv(path_out + "anno_cells_corr.txt"
                    ,sep = "\t"
                    ,index = True)



anno_samples=pd.concat([case, control])


# export anno_cells_corr
print("save anno_cells_corr.txt")
anno_samples.to_csv(path_out + "anno_samples_corr.txt"
                    ,sep = "\t"
                    ,index = True)
