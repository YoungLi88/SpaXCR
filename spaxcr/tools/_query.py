from anndata import AnnData
from scipy.sparse import csr_matrix
from itertools import chain
import numpy as np
from collections import defaultdict


from typing import Union, Literal, List, Optional


from .._utils import LiteralCheck, extract_var, GroupbyCheck


def query(
    adata: AnnData, 
    var_names: Optional[Union[str, List[str]]] = None, 
    obs_names: Optional[Union[str, List[str]]] = None, 
    obs_chain: Optional[str] = None,
    focus_cdr3: bool = False,
    cdr3_names: Optional[Union[str, List[str]]] = None, 
    cdr3_isotype: Optional[Literal['IGHA', 'IGHG', 'IGHM', 'IGHE', 'IGHD']] = None,
    flatten: bool = False,
    dedup: bool = True,
    remove_na: bool = False,
    reverse = False,
):
    """
    Query the AnnData object to return a dictionary of obs or var expressions based on the inputs.

    Args:
        adata (AnnData)
            AnnData object with a CSR matrix in adata.X.
        var_names
            One or more var names to query for their expressed obs.
        obs_names
            One or more obs names to query for their expressed vars.
        obs_chain:
            Only used for querying obs_names. Specific chains. Supports ['TRA','TRB','IGH','IGK','IGL'] or ['IGHA', 'IGHG', 'IGHM', 'IGHA', 'IGHD']
        focus_cdr3 : bool, optional, default=False
            Only used for querying obs_names. Whether to focus on CDR3 for expansion analysis. If True, the data will be grouped by CDR3.
        cdr3_names
            One or more obs names to query for their expressed cdr3.
        cdr3_isotype
            Only used for querying cdr3_names. Specific BCR isotype. Only supports ['IGHA', 'IGHG', 'IGHM', 'IGHA', 'IGHD']
        flatten
            Return List or Dict
    Returns:
        Dict: A dictionary where the keys are either obs names (if var_names is provided) 
              or var names (if obs_names is provided), and the values are lists of corresponding 
              expressed obs or var names.
    """
    # Ensure that the data matrix is CSR for efficiency
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)

    result = {}

    if var_names is not None:
        if isinstance(var_names, str):
            var_names = [var_names]
        
        for var in var_names:
            if var not in adata.var_names:
                raise ValueError(f"var '{var}' not found in adata.var_names")
            var_idx = adata.var_names.get_loc(var)
            # Find the obs indices where this var is expressed
            expressed_obs_indices = adata.X[:, var_idx].nonzero()[0]
            expressed_obs_names = adata.obs_names[expressed_obs_indices].tolist()
            result[var] = expressed_obs_names

    elif obs_names is not None:
        if isinstance(obs_names, str):
            obs_names = [obs_names]

        var_df = extract_var(adata, focus_cdr3 = focus_cdr3, chain = obs_chain)
        
        for obs in obs_names:
            obs_idx = adata.obs_names.get_loc(obs)
            # Find the var indices that are expressed in this obs
            expressed_var_indices = adata.X[obs_idx].nonzero()[1]
            if focus_cdr3:
                expressed_var_names = adata.var['CDR3'].iloc[expressed_var_indices].tolist()
            else:
                expressed_var_names = adata.var_names[expressed_var_indices].tolist()

            result[obs] = [ var for var in expressed_var_names if var in var_df['ID'].tolist() ]
            
    elif cdr3_names is not None:
        if isinstance(cdr3_names, str):
            cdr3_names = [cdr3_names]
            
        for cdr3 in cdr3_names:
            if cdr3 not in adata.var['CDR3'].tolist():
                raise ValueError(f"CDR3 '{cdr3}' not found in adata.var['CDR3']")
            
            if cdr3_isotype != None:
                valid_list = ['IGHA', 'IGHG', 'IGHM', 'IGHE', 'IGHD']
                LiteralCheck(cdr3_isotype, validlist = valid_list)
                found_loc = (adata.var['CDR3'] == cdr3) & (adata.var['IGHisotype'] == cdr3_isotype)
            else:
                found_loc = (adata.var['CDR3'] == cdr3)
            cdr3_idx = np.where(found_loc)[0]
            # Find the obs indices where this var is expressed
            expressed_obs_indices = adata.X[:,cdr3_idx].nonzero()[0]
            expressed_obs_names = adata.obs_names[expressed_obs_indices].tolist()
            result[cdr3] = expressed_obs_names

    if reverse and obs_names is None:
        reversed_dict = defaultdict(list)
        for key, values in result.items():
            for value in values:
                reversed_dict[value].append(key)
        result = reversed_dict

    if flatten :
        if dedup:
            result = list(set(chain.from_iterable(result.values())))
        else:
            result = list(chain.from_iterable(result.values()))
        return result
    
    if not remove_na:
        return result
    else:
        new = {}
        for k,v in result.items():
            if len(v) != 0:
                new[k] = v
        return new