from anndata import AnnData
import pandas as pd

from typing import Optional

from .._utils import n_cells_dict

def filter_clone(
    adata: AnnData,
    min_cdr3_length: int = 1,
    max_cdr3_length: int = 30,
    min_cells: int = 1,
    min_count: int = 1,
    count_column: str = None,
) -> AnnData:
    """
    Filters clones based on CDR3 length, number of cells, and count.

    Parameters:
        adata : AnnData
            Input AnnData object containing clone information.

        min_cdr3_length : int
            Minimum CDR3 length for filtering.

        max_cdr3_length : int
            Maximum CDR3 length for filtering.

        min_cells : int
            Minimum number of cells required for filtering.

        min_count : int
            Minimum count threshold for filtering clones.

        count_column : str
            The column name in clone_df representing the count to filter.


    Returns:
        AnnData: Filtered AnnData object after applying the filtering criteria.
    """

    # Copy the adata if replace is False
    adata = adata.copy()

    clone_df = adata.uns['immu_repertoire']['clone_info']
    clone_df['CDR3_length'] = clone_df['aaSeqCDR3'].map(len)
    raw_clone_num = len(adata.var_names)
    
    filter_clone_df = clone_df[
                (clone_df['CDR3_length'] >= min_cdr3_length) & 
                (clone_df['CDR3_length'] <= max_cdr3_length)
                ].copy()

    if min_count != 1:
        if count_column == None:
            candidate_col = [col for col in clone_df.columns if 'count' in col and 'mutated' not in col or 'Count' in col]
            if len(candidate_col) == 0:
                raise KeyError(f'Please specific columns of reads')
            count_column = candidate_col[0]
            print(f'Candidate columns of counts: {count_column}, Use {count_column}')

    # Check if the count_column exists, and filter clones based on criteria
        if count_column in filter_clone_df.columns:
            filter_clone_df = filter_clone_df[filter_clone_df[count_column] >= min_count].copy()
        else:
            candidate_col = [col for col in clone_df.columns if 'count' in col and 'mutated' not in col or 'Count' in col]
            raise KeyError(f'Not Found {count_column} in clone_df.columns. Candidate columns {candidate_col}')
    
    adata = adata[:, adata.var_names.isin(filter_clone_df['cloneId'])].copy()

    cells_dict = n_cells_dict(adata)

    adata.var['n_cells'] = adata.var_names.map(cells_dict)
    filtered_clone_num = len(adata.var_names)
#     print(f'Raw clones: {raw_clone_num}\nAfter filtered clones: {filtered_clone_num}')

    # filter cells (not suggest)
    adata = adata[:,adata.var['n_cells'] >= min_cells].copy()
    clone_df = clone_df[clone_df['cloneId'].isin(adata.var_names)].copy()
    adata.uns['immu_repertoire']['clone_info'] = clone_df

    return adata
