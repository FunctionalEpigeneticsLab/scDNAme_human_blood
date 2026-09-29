import os
import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import sys
sys.path.append('/data/leuven/357/vsc35706/notebook/sciMET/sciMET_PBMC/utils')
sys.path.append('/data/leuven/357/vsc35706/notebook/sciMET/scripts/analysis')
from utils import *
from conv_to_anndata import *

pkl_outdir = 'pkl'

def add_obs(adata):
    adata.obs['barcode'] = adata.obs.index.str.split('.var').str[0].str.split('.NCG').str[0]
    adata.obs.index.name = 'index'
    # Path to metadata
    meta_path = 'adata.singlets_only.processed.obs.csv'
    metadata = pd.read_csv(meta_path)
    metadata.index.name = 'index'
    print(metadata.shape)
    meta = adata.obs.merge(metadata, on='barcode', how='left')
    meta.index = meta['barcode'].copy()
    meta.index.name = 'index'
    adata.obs = meta.copy()
    adata_out = adata.copy()
    print(adata_out)
    return adata_out

# load original anndata object
input_data = {
              'PCA_Features_50k_Data_raw': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260120_BC_methscan/methscan_outs/50k_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasVMR_Data_raw': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260120_BC_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasCtgaVMR_Data_raw': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260120_BC_methscan/methscan_outs/VMR_debias_CT_GA_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasSnpVMR_Data_raw': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260120_BC_methscan/methscan_outs/VMR_commonSNP_matrix/matrices/MSR.h5ad',
              'PCA_Features_50k_Data_debias': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/50k_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasVMR_Data_debias': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasCtgaVMR_Data_debias': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/VMR_debias_CT_GA_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasSnpVMR_Data_debias': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/VMR_commonSNP_matrix/matrices/MSR.h5ad',
              'PCA_Features_50k_Data_debiasSNP': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_commonSNPs_methscan/methscan_outs/50k_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasVMR_Data_debiasSNP': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_commonSNPs_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasCtgaVMR_Data_debiasSNP': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_commonSNPs_methscan/methscan_outs/VMR_debias_CT_GA_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasSnpVMR_Data_debiasSNP': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_commonSNPs_methscan/methscan_outs/VMR_commonSNP_matrix/matrices/MSR.h5ad',
              'PCA_Features_50k_Data_debiasCtga': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_CT_GA_methscan/methscan_outs/50k_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasVMR_Data_debiasCtga': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_CT_GA_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasCtgaVMR_Data_debiasCtga': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_CT_GA_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad',
              'PCA_Features_debiasSnpVMR_Data_debiasCtga': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_CT_GA_methscan/methscan_outs/VMR_commonSNP_matrix/matrices/MSR.h5ad',
              'PCA_Features_rawVMR_Data_raw': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260120_BC_methscan/methscan_outs/VMR_raw_matrix/matrices/MSR.h5ad',
              'PCA_Features_rawVMR_Data_debias': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/VMR_raw_matrix/matrices/MSR.h5ad',
              'PCA_Features_rawVMR_Data_debiasCtga': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_CT_GA_methscan/methscan_outs/VMR_raw_matrix/matrices/MSR.h5ad',
              'PCA_Features_rawVMR_Data_debiasSNP': '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_commonSNPs_methscan/methscan_outs/VMR_raw_matrix/matrices/MSR.h5ad',
              }

# PIPELINE to prepare data for scIB
for prefix, adata_file in input_data.items():
    adata = sc.read_h5ad(adata_file)

    add_obs(adata)

    mask = ((adata.obs["cell_type_main"].notna()))
    adata = adata[mask].copy()
    print(adata)

    sc.tl.pca(adata, zero_center=True)

    obsm_obj = adata.obsm["X_pca"].copy()[:, :15]
    obsm_di = {prefix: obsm_obj}
    outfile = os.path.join(pkl_outdir, f'{prefix}.obsm_di.pkl')
    with open(outfile, "wb") as f:
        pickle.dump(obsm_di, f)

    n_components = 15
    n_neighbors = 20
    print(f'Using {n_components} PCs to compute NN graph.')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_components)
    sc.tl.umap(adata, random_state=42, min_dist=0.15)

    # save figures
    plt.rcParams['figure.figsize'] = (3.3,3)
    alpha = 1

    sc.pl.umap(adata, color='cell_type_main', title = 'Cell type', palette = celltype_main_palette,
               alpha=alpha, legend_loc='on data', legend_fontsize='xx-small', legend_fontweight='normal', frameon=False,
               save=f'.{prefix}.color_by_cell_type_main.pdf')

    sc.pl.umap(adata, color='cell_type_main', title = 'Cell type', palette = celltype_main_palette,
               alpha=alpha, frameon=False,
               save=f'.{prefix}.color_by_cell_type_main.sep_legend.pdf')

    sc.pl.umap(adata, color='donor_id', title='Donor', palette = donor_palette2,
               alpha=alpha, frameon=False,
               save=f'.{prefix}.color_by_donor_id.pdf')

    sc.pl.umap(adata, color='cell_type_main', title = 'Cell type', palette = celltype_main_palette,
               alpha=alpha, frameon=False,
               save=f'.{prefix}.color_by_cell_type_main.sep_legend.png')

    sc.pl.umap(adata, color='donor_id', title='Donor', palette = donor_palette2,
               alpha=alpha, frameon=False,
               save=f'.{prefix}.color_by_donor_id.png')

    sc.pl.umap(adata, color='experiment', alpha=alpha, frameon=False, palette = experiment_palette,
               save=f'.{prefix}.color_by_experiment.pdf')

    sc.pl.umap(adata, color='NCGN_perc', title='%mCG', alpha=alpha, frameon=False, save=f'.{prefix}.color_by_mCG.png')
    sc.pl.umap(adata, color='NCHN_perc', title='%mCH', alpha=alpha, frameon=False, save=f'.{prefix}.color_by_mCH.png')
    sc.pl.umap(adata, color='coverage', title='% genome-wide coverage', alpha=alpha, vmax=5, frameon=False, save=f'.{prefix}.color_by_coverage.png')

# Save a final anndata object for scIB input
adata_file = '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/methscan_outs/VMR_debias_matrix/matrices/MSR.h5ad'
adata = sc.read_h5ad(adata_file)
adata.obsm.pop("X_pca", None)
adata.layers.pop('MSR', None)

add_obs(adata)

# Filter object
mask = ((adata.obs["cell_type_main"].notna()))
adata = adata[mask].copy()
print(adata)

# Load computed PCA obsm
pkl_files = glob.glob(os.path.join(pkl_outdir, "*.obsm_di.pkl"))
print('Number of pkl files:')
print(len(pkl_files))

for pkl_file in pkl_files:
    with open(pkl_file, "rb") as f:
        obsm_di = pickle.load(f)
        for prefix, obsm_obj in obsm_di.items():
            adata.obsm[prefix] = obsm_obj.copy()
print(adata.obsm)

# save data
def has_bad_objects(s):
    return s.dtype == object and not s.dropna().map(lambda x: isinstance(x, str)).all()

for df in (adata.obs, adata.var):
    bad = [c for c in df.columns if has_bad_objects(df[c])]
    print('dropping:', bad)
    df.drop(columns=bad, inplace=True)

ad.settings.allow_write_nullable_strings = True
adata_file = 'adata/adata_for_scIB.h5ad'
adata.write(adata_file)
