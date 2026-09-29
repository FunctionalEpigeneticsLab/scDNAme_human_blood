import argparse
import mudata
import numpy as np
import scanpy as sc
import scvi
import torch
from scvi.external import METHYLVI
import anndata as ad

def main():
    parser = argparse.ArgumentParser(description='Run METHYLVI analysis on mudata object')
    parser.add_argument('--input', required=True, help='Input mudata file (.h5mu)')
    parser.add_argument('--output', required=True, help='Output mudata file (.h5mu)')
    parser.add_argument('--n_latent', type=int, default=20, help='Number of latent dimensions')
    parser.add_argument('--max_epochs', type=int, default=500, help='Maximum number of training epochs')
    parser.add_argument('--early_stopping', action='store_true', default=True, help='Enable early stopping (default: True)')
    parser.add_argument('--no-early_stopping', dest='early_stopping', action='store_false', help='Disable early stopping')
    parser.add_argument('--batch_key', default='experiment')
    args = parser.parse_args()
    mdata = mudata.read(args.input)
    print(f"Shape: {mdata.shape}")
    
    METHYLVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        methylation_contexts=["mCG"],
        categorical_covariate_keys=[args.batch_key],
        modalities={"categorical_covariate_keys": "mCG"},
    )
    
    model = METHYLVI(
        mdata,
        n_latent=args.n_latent,
    )
    
    model.train(max_epochs=args.max_epochs, early_stopping=args.early_stopping)
    
    mdata.obsm["methylVI"] = model.get_latent_representation()
    
    sc.pp.neighbors(mdata, use_rep="methylVI")
    
    ad.settings.allow_write_nullable_strings = True
    mdata.write(args.output)

if __name__ == '__main__':
    main()
