import os
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import sys
sys.path.append('/data/leuven/357/vsc35706/notebook/sciMET/sciMET_PBMC/utils')
sys.path.append('/data/leuven/357/vsc35706/notebook/sciMET/scripts/analysis')
from conv_to_anndata import *

adata_outdir = '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_snm3Cseq_integration/adata'
os.makedirs(adata_outdir, exist_ok=True)

# get metadata
adata_file = '/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan/adata/adata.singlets_only.processed.h5ad'
adata = sc.read_h5ad(adata_file)
adata.obs.drop(columns=['keep'], inplace=True)
metadata = adata.obs.copy()
del(adata)

input_dir = os.path.join('/staging/leuven/stg_00064/Xinran/project/sciMET_PBMC/260918_BC_debias_methscan', 'methscan_outs', 'VMR_matrix/matrices')
mdata = methscan_outs_to_mudata(input_dir)
del mdata['mCG'].layers['MSR']

# Add metadata
def add_obs(adata, metadata):
    metadata.index.name = 'index'
    adata.obs['barcode'] = adata.obs.index.str.split('.var').str[0].str.split('.NCG').str[0]
    adata.obs.index.name = 'index'
    meta = adata.obs.merge(metadata, on='barcode', how='left')
    meta.index = meta['barcode'].copy()
    meta.index.name = 'index'
    adata.obs = meta.copy()
    adata_out = adata.copy()
    return adata_out

add_obs(mdata['mCG'], metadata)

print(mdata['mCG'].obs)

# Filter object
mask = (
    (mdata['mCG'].obs["cell_type_main"].notna())
)

mdata = mdata[mask].copy()
mdata

# save mudata to disk
import anndata as ad
ad.settings.allow_write_nullable_strings = True
prefix = 'sciMET_singlets_only'
outfile = os.path.join(adata_outdir, f"mdata_input.{prefix}.h5mu")
mdata.write(outfile)

mdata_inhouse = mdata.copy()
del(mdata)
mdata_inhouse.obs['DATASET'] = 'sciMET'
mdata_inhouse['mCG'].obs['DATASET'] = 'sciMET'

###################

# load ref data
mdata = mudata.read('/staging/leuven/stg_00064/Xinran/project/snm3Cseq/260130_HumanCellEpigenomeAtlas_PBMC_6082cells_methscan/adata/mdata_input.VMR_debias.h5mu')
input_dir = os.path.join('/staging/leuven/stg_00064/Xinran/project/snm3Cseq/260130_HumanCellEpigenomeAtlas_PBMC_6082cells_methscan', 'methscan_outs', 'VMR_debias_matrix/matrices')
mdata = methscan_outs_to_mudata(input_dir)
del mdata['mCG'].layers['MSR']

# process obs
metadata = mdata['mCG'].obs.copy()
metadata['filename'] = metadata.index
metadata.index.name = 'index'
metadata['barcode'] = metadata['filename'].str.split('.CGN').str[0]
# add batch key
metadata['batch'] = metadata['barcode'].str.split('_Plate').str[0].astype("category")
metadata['DATASET'] = 'snm3Cseq'
mdata['mCG'].obs = metadata.copy()
mdata['mCG'].var.index.name = 'index'
mdata.obs.index.name = 'index'

# load metadata
meta_file = '/staging/leuven/stg_00064/publicdata/20261201_HumanCellEpigenomeAtlas_sc_allc_Ecker/metadata/HumanCellEpigenomeAtlas_metadata_with_anno.csv'
meta_all = pd.read_csv(meta_file)
meta = meta_all.loc[meta_all['tissue_type']=='PBMC']
meta.drop(columns = ['Cis/Trans', 'Short/Long'], inplace=True)
meta.shape[0]
# update metadata
meta_updated = metadata.merge(meta, on='barcode', how='left')
meta_updated.index = meta_updated['barcode']
meta_updated.index.name = 'index'
mdata['mCG'].obs = meta_updated.copy()
mdata.obs = meta_updated.copy()

# filter
excl_celltypes = [
    "Hema NK CD16 Lung",
    "Hema Myeloid Monocyte Heart",
    "Hema Tmem Effector CD8 Mix",
    "Hema B Plasma",
    "Hema NK CD56 Mix",
    "Hema B Memory 1",
    "Hema NK ILC",
    "Hema Tmem Regulatory CD4 Skin",
    "Hema NK CD56 Skin",
    "Hema Tmem Follicular Helper CD4 Intestinal",
    "Hema Tmem Gamma Delta CD8 Lung",
]

mask = (
    mdata.obs["celltype_L2_both_annot"].notna()
    & ~mdata.obs["celltype_L2_both_annot"].isin(excl_celltypes)
)

mdata = mdata[mask].copy()
mdata

# save mudata to disk
import anndata as ad
ad.settings.allow_write_nullable_strings = True
prefix = 'snm3Cseq'
outfile = os.path.join(adata_outdir, f"mdata_input.{prefix}.h5mu")
mdata.write(outfile)


#######################

# subset mdata_inhouse to overlapping vars
ref_vars = mdata["mCG"].var_names
adata_inhouse = mdata_inhouse["mCG"]
common = adata_inhouse.var_names.isin(ref_vars)
adata_inhouse._inplace_subset_var(common)

# make a new mdata out of anndata
mdata_inhouse = mudata.MuData({"mCG": adata_inhouse})

# Concat two datasets
mdata_concat = mudata.concat(
    [mdata, mdata_inhouse],
    join="outer"
)
mdata_concat

# set batch key
del mdata_concat["mCG"].obsm["X_pca"]
obs = mdata_concat["mCG"].obs
obs["experiment"] = obs["experiment"].astype(str)
obs["donor_id"] = obs["donor_id"].astype(str)
obs["batch"] = obs["batch"].astype(str)
obs["experiment"] = obs["experiment"].fillna(obs["batch"])
print(obs["experiment"].value_counts())

donor = obs["donor_id"].astype(str)
batch = obs["batch"].astype(str)
mask = donor.isin(["nan", "NA", "None", ""]) | obs["donor_id"].isna()
obs["donor_id"] = donor.where(~mask, batch)
obs["donor_id"] = obs["donor_id"].astype("category")
mdata.update()
print(obs["donor_id"].value_counts())

# modify
mdata_concat.obs = obs.copy()
mdata_concat['mCG'].obs = obs.copy()

# save mdata
mdata_file = os.path.join(adata_outdir, 'mdata_input.VMR_debias.concat_snm3Cseq_sciMET.h5mu')
# save
ad.settings.allow_write_nullable_strings = True
mdata_concat.write(mdata_file)

