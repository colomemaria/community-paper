#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scanpy as sc
import pandas as pd
import os
import shutil
import sys

def process_data(n_case, n_control, path_in, path_out):
    data = sc.read_csv(path_in + "counts_corr.csv.gz")
    anno_cells = pd.read_csv(path_in + "anno_cells_corr.txt", sep="\t")
    anno_samples = pd.read_csv(path_in + "anno_samples_corr.txt", sep="\t")

    case = anno_samples[anno_samples.case_or_control == "case"]
    control = anno_samples[anno_samples.case_or_control == "control"]
    data = data.T
    anno_cells.index = anno_cells.cell_ID
    data.obs = anno_cells.loc[data.obs.index]

    case = case.sample(n_case)
    control = control.sample(n_control)

    ID_list = control.sample_ID.to_list() + case.sample_ID.to_list()
    sub = data[data.obs.sample_ID.isin(ID_list)].copy()
    sub = sub.T

    output_dir = os.path.join(path_out, f"{n_case}_{n_control}")
    os.makedirs(output_dir, exist_ok=True)

    sub.to_df().to_csv(os.path.join(output_dir, "counts_corr.csv.gz"), index=True, compression="gzip")

    # export anno_cells_corr
    print("save anno_cells_corr.txt")
    sub.var.to_csv(os.path.join(output_dir, "anno_cells_corr.txt"), sep="\t", index=True)

    anno_samples = pd.concat([case, control])

    # export anno_samples_corr
    print("save anno_samples_corr.txt")
    anno_samples.to_csv(os.path.join(output_dir, "anno_samples_corr.txt"), sep="\t", index=True)

if __name__ == "__main__":
    path_in = str(sys.argv[1])
    path_out = str(sys.argv[2])

    case_values = [3, 7, 10, 17]
    control_values = [3, 6, 12, 11]

    for n_case, n_control in zip(case_values, control_values):
        process_data(n_case, n_control, path_in, path_out)

