import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import sys
sys.path.append('/data/leuven/357/vsc35706/notebook/sciMET/scripts/analysis')
from doublet_detection import *

adata_file = 'adata/adata.processed.h5ad'
adata = sc.read_h5ad(adata_file)

mask = adata.obs["donor_id"] != "doublet"
adata_sub = adata[mask].copy()
print(adata_sub)
del(adata)

doublet_adata = create_artificial_doublets(adata_sub, sim_fraction=1, matrix_keys=['MSR'], return_simulated_only=True, n_jobs=16)
print(doublet_adata)

adata_file = 'adata/adata.processed.VMR.artificial_doublets.1.h5ad'
doublet_adata.write(adata_file)
