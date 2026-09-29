import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from joblib import Parallel, delayed
from typing import Optional, List
from sklearn.neighbors import NearestNeighbors

def create_artificial_doublets(
    adata: anndata.AnnData,
    sim_fraction: float = 1,
    matrix_keys: Optional[List[str]] = None,
    n_jobs: int = -1,
    random_state: int = 42,
    return_simulated_only: bool = False
) -> anndata.AnnData:
    np.random.seed(random_state)

    n_cells = adata.n_obs
    n_doublets = int(np.ceil(n_cells * sim_fraction))
    cell_names = np.array(adata.obs_names)
    pairing_names = np.random.choice(cell_names, size=(2, n_doublets), replace=True)

    if matrix_keys is None:
        matrix_keys = list(adata.layers.keys()) if len(adata.layers) > 0 else ['X']

    obs_names_combined = list(adata.obs_names) + [f'Doublet_{i}' for i in range(n_doublets)]
    obs_names_sim = [f'Doublet_{i}' for i in range(n_doublets)]

    new_obs_combined = pd.DataFrame(
        {'doublet_status': ['True_Cell'] * n_cells + ['Artificial_Doublet'] * n_doublets},
        index=obs_names_combined
    )
    new_obs_sim = pd.DataFrame(
        {'doublet_status': ['Artificial_Doublet'] * n_doublets},
        index=obs_names_sim
    )

    result_dict = {} 

    for matrix_key in matrix_keys:
        if matrix_key == 'X':
            matrix = adata.X.copy()
        else:
            matrix = adata.layers[matrix_key].copy()
        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()
        if matrix.shape[0] != adata.n_vars:
            matrix = matrix.T

        def create_doublet_column(col_idx):
            cell1_name = pairing_names[0, col_idx]
            cell2_name = pairing_names[1, col_idx]
            cell1_idx = adata.obs_names.get_loc(cell1_name)
            cell2_idx = adata.obs_names.get_loc(cell2_name)
            vals = np.array([matrix[:, cell1_idx], matrix[:, cell2_idx]])
            result = np.zeros(adata.n_vars)
            for i in range(adata.n_vars):
                cell_vals = vals[:, i]
                if np.all(np.isnan(cell_vals)):
                    result[i] = np.nan
                elif np.any(np.isnan(cell_vals)):
                    if np.random.random() < 0.5:
                        result[i] = np.nanmean(cell_vals)
                    else:
                        result[i] = np.nan
                else:
                    result[i] = np.mean(cell_vals)
            return result
        doublet_mat = np.column_stack(
            Parallel(n_jobs=n_jobs)(delayed(create_doublet_column)(i) for i in range(n_doublets))
        )
        if return_simulated_only:
            result_dict[matrix_key] = doublet_mat.T
        else:
            combined_mat = np.hstack([matrix, doublet_mat])
            result_dict[matrix_key] = combined_mat.T 

    import anndata
    if return_simulated_only:
        main_key = matrix_keys[0]
        doublet_adata = anndata.AnnData(
            X=result_dict[main_key],
            obs=new_obs_sim,
            var=adata.var.copy()
        )
        for k in matrix_keys:
            if k != 'X':
                doublet_adata.layers[k] = result_dict[k]
    else:
        main_key = matrix_keys[0]
        doublet_adata = anndata.AnnData(
            X=result_dict[main_key],
            obs=new_obs_combined,
            var=adata.var.copy()
        )
        for k in matrix_keys:
            if k != 'X':
                doublet_adata.layers[k] = result_dict[k]
        # Add metadata from original (same as your logic)
        for col in adata.obs.columns:
            if col not in doublet_adata.obs.columns:
                doublet_adata.obs[col] = adata.obs[col].tolist() + [np.nan]*n_doublets

    return doublet_adata


def compute_doublet_scores(
    X_cells, 
    X_doublets,
    n_neighbors=None, 
    random_state=0
):
    n_cells = X_cells.shape[0]
    n_doublets = X_doublets.shape[0]
    total = n_cells + n_doublets

    if n_neighbors is None:
        n_neighbors = int(round(0.5 * np.sqrt(n_cells)))

    X_all = np.vstack([X_cells, X_doublets])
    labels = np.concatenate([
        np.zeros(n_cells, dtype=bool),
        np.ones(n_doublets, dtype=bool),
    ])

    nn = NearestNeighbors(n_neighbors=n_neighbors + 1,
                          algorithm='auto')
    nn.fit(X_all)
    dists, nbrs = nn.kneighbors(X_all)
    is_doublet = labels[nbrs] 
    is_doublet = is_doublet[:, 1:]
    doublet_neighbor_count = is_doublet.sum(axis=1)
    fraction_doublet_neighbors = doublet_neighbor_count / n_neighbors
    doublet_scores_cells = fraction_doublet_neighbors[:n_cells]
    doublet_scores_doublets = fraction_doublet_neighbors[n_cells:]
    return doublet_scores_cells, doublet_scores_doublets



def compute_doublet_scores_3group(
    X_real,  
    X_doublets, 
    X_heterogenic,
    n_neighbors=None,
    random_state=0
):

    n_real = X_real.shape[0]
    n_doublets = X_doublets.shape[0]
    n_het = X_heterogenic.shape[0]

    if n_neighbors is None:
        n_neighbors = int(round(0.5 * np.sqrt(n_real)))

    X_all = np.vstack([X_real, X_doublets, X_heterogenic])

    labels_is_doublet = np.concatenate([
        np.zeros(n_real, dtype=bool), 
        np.ones(n_doublets, dtype=bool),
        np.zeros(n_het, dtype=bool),
    ])

    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, algorithm='auto')
    nn.fit(X_all)
    dists, nbrs = nn.kneighbors(X_all)
    is_doublet_neighbor = labels_is_doublet[nbrs]
    is_doublet_neighbor = is_doublet_neighbor[:, 1:]  
    doublet_neighbor_count = is_doublet_neighbor.sum(axis=1)
    fraction_doublet_neighbors = doublet_neighbor_count / n_neighbors
    doublet_scores_real = fraction_doublet_neighbors[:n_real]
    doublet_scores_doublets = fraction_doublet_neighbors[n_real:n_real + n_doublets]
    doublet_scores_heterogenic = fraction_doublet_neighbors[n_real + n_doublets:]
    return doublet_scores_real, doublet_scores_doublets, doublet_scores_heterogenic
