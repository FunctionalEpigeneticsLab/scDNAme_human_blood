import pandas as pd
import numpy as np
import scipy.sparse as sps
import scipy.io as io
import scanpy as sc
import gzip
import os
import argparse

def csv_to_sp(input_file):

    def _process_data(row_array):
        return row_array
    
    sp_data = []
    row_names = []
    col_names = None
    sep = ','
    
    with gzip.open(input_file, 'rt') as csv_file:
        for i, row in enumerate(csv_file):
            if i == 0:
                col_names = row.strip().split(sep)[1:] 
                continue
            
            # parse one row at a time
            cols = row.strip().split(sep)
            row_name = cols[0]
            # convert empty strings to nan 
            data = np.array([float(x) if x != '' else np.nan for x in cols[1:]])
            
            row_names.append(row_name)
            data = _process_data(data)
            data = sps.csr_matrix(data)
            sp_data.append(data)
    
    sp_data = sps.vstack(sp_data).tocsr()

    # remove NaN values from sparse matrix
    sp_data_clean = sp_data.copy()
    sp_data_clean.data[np.isnan(sp_data_clean.data)] = 0
    sp_data_clean.eliminate_zeros()

    print(f"Shape of matrix: {sp_data_clean.shape}")
    print(f"Row names: {row_names[:5]}")
    print(f"Column names:  {col_names[:5]}")
    return sp_data_clean, row_names, col_names

def write_mm(sp_data, row_names, col_names, outdir='matrices', prefix='MSR', write_cells_and_features=True, write_anndata=True):

    os.makedirs(outdir, exist_ok=True)

    # write Matrix Market format (. mtx)
    import scipy.io as io
    io.mmwrite(os.path.join(outdir, f'{prefix}.mtx'), sp_data)  # cells x features
    print("mtx file created.")
    
    if write_cells_and_features:
        barcodes_df = pd.DataFrame({'cell_id': row_names})
        barcodes_df.to_csv(os.path.join(outdir, 'barcodes.tsv'), sep='\t', header=False, index=False)
        
        features_df = pd.DataFrame({
            'feature_id': col_names
        })
        features_df.to_csv(os.path.join(outdir, 'features.tsv'), sep='\t', header=False, index=False)
        print("Barcodes and features TSV files created.")
    
    if write_anndata: 
        import anndata as ad
        # create AnnData object (cells x features)
        adata = ad.AnnData(
            X=sp_data,
            obs=pd.DataFrame(index=row_names),
            var=pd.DataFrame(index=col_names)
        )
        
        # freeze layer
        adata.layers['MSR'] = adata.X.copy()
        
        # truncated SVD
        sc.tl.pca(adata, zero_center=True)
        
        adata.write_h5ad(os.path.join(outdir, f'{prefix}.h5ad'))
        print('AnnData object created.')


def main():
    parser = argparse.ArgumentParser(description='Convert csv.gz (MethSCAn output) to Matrix Market and/or AnnData format')
    
    parser.add_argument('--input_file', required=True, help='Input csv.gz')
    parser.add_argument('--outdir', default='matrices', help='Output directory (default: matrices)')
    parser.add_argument('--prefix', default='MSR', help='Output file prefix (default: MSR)')
    parser.add_argument('--write_cells_and_features', action='store_true', default=False,
                        help='Write barcodes and features tsv (default: False)')
    parser.add_argument('--write_anndata', action='store_true', default=False,
                        help='Write AnnData h5ad file (default: False)')
    
    args = parser.parse_args()
    
    print(f"Processing:  {args.input_file}")
    sp_data, row_names, col_names = csv_to_sp(args.input_file)
    
    write_mm(sp_data, row_names, col_names, 
             outdir=args.outdir, 
             prefix=args.prefix,
             write_cells_and_features=args.write_cells_and_features,
             write_anndata=args.write_anndata)

if __name__ == '__main__':
    main()
