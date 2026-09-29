import os
import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pickle

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from dataclasses import fields
print([f.name for f in fields(BatchCorrection)])
print([f.name for f in fields(BioConservation)])


# read saved anndata object
adata_file = 'adata/adata_for_scIB.h5ad'
adata = sc.read_h5ad(adata_file)

from itertools import product

feature_order = [
    "rawVMR",
    "debiasCtgaVMR",
    "debiasVMR",
    "debiasSnpVMR",
    "50k",
]

data_order = [
    "raw",
    "debiasCtga",
    "debias",
    "debiasSNP",
]

# Features first, then Data within each Feature
embedding_obsm_keys = [
    f"PCA_Features_{feature}_Data_{data}"
    for feature in feature_order
    for data in data_order
]

for i, key in enumerate(embedding_obsm_keys, start=1):
    print(f"emb{i:02d}: {key}")

embedding_map = {
    f"emb{i:02d}": key
    for i, key in enumerate(embedding_obsm_keys, start=1)
}

for emb_name, original_key in embedding_map.items():
    adata.obsm[emb_name] = adata.obsm[original_key]


for emb_name, original_key in embedding_map.items():
    print(f"{emb_name} = {original_key}")

bm = Benchmarker(
    adata,
    batch_key="donor_id",
    label_key="cell_type_main",
    bio_conservation_metrics=BioConservation(
        clisi_knn=False
    ),
    batch_correction_metrics=BatchCorrection(
        bras=False,
        pcr_comparison=False,
    ),
    embedding_obsm_keys=list(embedding_map.keys()),
    n_jobs=32,
)


bm.benchmark()

# save results
outfile = 'pkl/bm.pkl'
with open(outfile, "wb") as f:
        pickle.dump(bm, f)

# save figure
bm.plot_results_table(min_max_scale=False, save_dir="figures")

from rich import print
df = bm.get_results(min_max_scale=False)
print(df)
df.to_csv('bm.results.final.csv')

