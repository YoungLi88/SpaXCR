from anndata import AnnData
import pandas as pd
import numpy as np
import scanpy as sc

from typing import Optional


def addloc(
    adata: AnnData,
    loc_info: pd.DataFrame,
    basis: str = 'spatial',
) -> AnnData:
    """
    Adds x, y location information to the AnnData object.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object to which the location information will be added. 
        
    loc_info : pandas.DataFrame
        A DataFrame containing location information with barcodes as the index.
        The DataFrame must include 'x' and 'y' columns representing spatial coordinates.


    basis: str
        obsm contains Location information
    
    Returns:
    --------
    AnnData
        The generated AnnData object containing the location information.
    """

    adata = adata.copy()

    if not isinstance(loc_info, pd.DataFrame):
        raise ValueError("loc_info must be a pandas DataFrame.")
    
    required_columns = {'x', 'y'}
    if not required_columns.issubset(loc_info.columns):
        raise ValueError("loc_info DataFrame must contain 'x' and 'y' columns.")
    
    common_index = adata.obs_names.intersection(loc_info.index)

    if len(common_index) == 0:
        raise ValueError('Not Found any intersection between adata.obs_names and loc_info.index')

    adata = adata[common_index]
    rest_cells = set(loc_info.index) - set(adata.obs_names)
    rest_df = pd.DataFrame(
        np.zeros((len(rest_cells),adata.shape[1])),
        index = list(rest_cells),
        columns = adata.var_names
        )
    
    rest_adata = AnnData(rest_df)

    final_adata = sc.concat([adata,rest_adata],axis = 0)
    final_adata = final_adata[loc_info.index].copy()
    final_adata.var = adata.var
    final_adata.obs['x'] = loc_info['x']
    final_adata.obs['y'] = loc_info['y']
    final_adata.obsm[basis] = final_adata.obs[['x','y']].values
    if 'spatial' in adata.uns.keys():
        final_adata.uns['spatial'] = adata.uns['spatial']
    final_adata.uns['immu_repertoire'] = adata.uns['immu_repertoire']
    final_adata = final_adata.copy()
    return final_adata


