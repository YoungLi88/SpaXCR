from typing import Union, List, Optional
import pandas as pd
import numpy as np
from anndata import AnnData

from sklearn.neighbors import KDTree

from .._utils import LiteralCheck, n_cells_dict, extract_var



def expansion(
    adata: AnnData,
    chain: Optional[str] = None,
    focus_cdr3: bool = False,
    groups: Union[List[int], int] = 1,
    labels: Optional[List[str]] = ['Notexpansion', 'Expansion'],
    show_match_dict: bool = True,
) -> pd.DataFrame:
    """
    Categorize cells based on the number of cells per CDR3 or chain.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing clone information.
    
    chain : str
        Chain type to focus on. Should be one of ['TRA','TRB','IGH','IGK','IGL','IGHA','IGHG','IGHM','IGHD','IGHE'].
    
    focus_cdr3 : bool, optional, default=False
        Whether to focus on CDR3 for expansion analysis. If True, the data will be grouped by CDR3.

    groups : Union[List[int], int], optional, default=1
        Defines the intervals of cell counts to categorize CDR3 or chains. 
        If groups is an integer, divides into (0-groups, >groups).
        If groups is a list, divides into intervals (0, list[0], list[1], ..., >list[-1]).

    labels : List[str], optional, default=['Notexpansion', 'Expansion']
        List of labels for the categories defined by `groups`. The number of labels must match the number of intervals.

    show_match_dict : bool, optional, default = False
        Show match dict or not

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the sorted and categorized information for CDR3/chain with the associated expansion groups.
    """


    var_df = extract_var(adata, chain = chain, focus_cdr3 = focus_cdr3)


    if chain == None or len(chain) == 3 :
        sort_chainstat = var_df.sort_values('n_cells', ascending=False, ignore_index=True)
        sort_chainstat = sort_chainstat[['ID', 'Chain', 'n_cells']].copy()
    else:
        sort_chainstat = var_df.sort_values('n_cells', ascending=False, ignore_index=True)
        sort_chainstat = sort_chainstat[['ID', 'IGHisotype', 'n_cells']].copy()


    sort_chainstat, group_labels = categorize_n_cells(sort_chainstat, groups=groups)

    # If labels are provided, ensure they match the number of groups and apply them
    if labels is not None:
        if len(labels) != len(group_labels):
            raise KeyError(f'labels do not match group labels. Groups contain {group_labels}.')
        else:
            match_dict = dict(zip(group_labels, labels))
            if show_match_dict:
                print(f'Match Dict: {match_dict}')
            sort_chainstat['group'] = sort_chainstat['group'].map(match_dict)

    return sort_chainstat

def categorize_n_cells(
    df: pd.DataFrame,
    groups: Union[List,int] = 1,
):
    if isinstance(groups,int) or len(groups) == 1:
        if isinstance(groups,int):
            groups = [groups]
        bins = [-float('inf')] + groups+ [float('inf')]
        if groups[0] == 1:
            start = '1'
        else:
            start = f'<={groups[0]}'
        labels = [start] + [f'>{groups[0]}'] 
        
        
    else:
        groups  = sorted((groups))
        bins = [-float('inf')] + groups + [float('inf')]
        if groups[0] == 1:
            start = '1'
        else:
            start = f'<={groups[0]}'
        labels = [start] + [f'{bins[i]}-{bins[i+1]}' for i in range(1, len(bins) - 2)] \
        + [f'>{bins[-2]}']
    
    df['group'] = pd.cut(df['n_cells'], bins = bins, labels = labels, right = True)
    return df, labels
