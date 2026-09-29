import glob
import os
import scanpy as sc
from scipy import io
import mudata
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

def anndata_to_mudata(adata, name):
    # create mudata object
    mdata = mudata.MuData({name: adata})
    mdata.var.index.name='index'
    return mdata
    
def methscan_outs_to_mudata(input_dir, name="mCG"):
    mtx_cov = io.mmread(os.path.join(input_dir, "cov.mtx")).tocsr()
    mtx_cov.data = mtx_cov.data.astype(int)
    mtx_mc = io.mmread(os.path.join(input_dir, 'mc.mtx')).tocsr()
    mtx_mc.data = mtx_mc.data.astype(int)
    adata_file = os.path.join(input_dir,'MSR.h5ad')
    adata = sc.read_h5ad(adata_file)
    adata.layers['cov'] = mtx_cov
    adata.layers['mc'] = mtx_mc
    
    mdata = anndata_to_mudata(adata, name)
    return mdata

def adata_to_df(adata, cov_cutoff=10):
    frac = adata.layers['frac']
    cov = adata.layers['cov']

    frac = frac.toarray() if sparse.issparse(frac) else np.asarray(frac)
    cov = cov.toarray() if sparse.issparse(cov) else np.asarray(cov)

    df = pd.DataFrame(frac, index=adata.obs.index, columns=adata.var.index)
    cov_df = pd.DataFrame(cov, index=adata.obs.index, columns=adata.var.index)

    df = df.where(cov_df >= cov_cutoff)

    return df

def clean_up_h5ad(adata_file):
    import h5py
    
    problematic_keys = ['log1p', 'dendrogram_leiden', 'neighbors', 'rank_genes_groups', 
                         'hvg', 'pca', 'umap']
    
    with h5py.File(adata_file, 'r+') as f:
        if 'uns' in f:
            for key in problematic_keys:
                if key in f['uns']:
                    try:
                        del f[f'uns/{key}']
                        print(f"Removed /uns/{key}")
                    except Exception as e:
                        print(f"Failed to remove /uns/{key}: {e}")
