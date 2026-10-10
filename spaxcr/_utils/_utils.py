from anndata import AnnData
from typing import Optional, List


def LiteralCheck(inele,validlist,valid_name = None,return_detail = True):

    if isinstance(inele,str):
        if inele not in validlist:
            raise KeyError(f'Invalid {inele} keys. Only supports {validlist}')
    elif isinstance(inele,list):
        invalid_list = []
        for ele in inele:
            if ele not in validlist:
                invalid_list.append(ele)
        if len(invalid_list) != 0:
            if return_detail:
                raise KeyError(f'Invalid {invalid_list} keys. Only supports {validlist}')
            else:
                raise KeyError(f'Invalid {invalid_list} keys in {valid_name}')

def valid_chain(
    adata: AnnData,
    isotype: bool = False,
):
    chains = adata.var['Chain'].unique().tolist()
    if isotype and 'IGHisotype' in adata.var.columns:
        isochain = adata.var['IGHisotype'].value_counts().index.tolist()
        chains += isochain
    return chains

def n_cells_dict(
    adata: AnnData
):
    cells_dict = {}
    for i, var_name in enumerate(adata.var_names):
        expressed_obs_indices = adata.X[:, i].nonzero()[0]
        cells_num = len(expressed_obs_indices)
        cells_dict[var_name] = cells_num
    return cells_dict


def extract_var(
    adata: AnnData,
    focus_cdr3: bool = False,
    chain: Optional[str] = None,
):

    cells_dict = n_cells_dict(adata)
    adata.var['n_cells'] = adata.var_names.map(cells_dict)
    var_df = adata.var.copy()

    # subset clone_type
    if chain != None:
        
        if len(chain) == 3:
            
            valid_list = var_df['Chain'].unique().tolist()
            LiteralCheck(chain, valid_list)
            var_df = var_df[var_df['Chain'] == chain].copy()

        elif len(chain) == 4:
            if 'IGHisotype' not in var_df.columns:
                raise KeyError('Not Found IGHisotype in adata.var. Please enrich your clone information')

            valid_list = var_df['IGHisotype'].unique().tolist()
            LiteralCheck(chain, valid_list)
            var_df = var_df[var_df['IGHisotype'] == chain].copy()
        else:
            raise KeyError(f'Invalid chain {chain}')

    var_df['CDR3'] = var_df['CDR3'].astype(object)
    var_df['Chain'] = var_df['Chain'].astype(object)
    if 'IGHisotype' in var_df.columns:
        var_df['IGHisotype'] = var_df['IGHisotype'].astype(object).fillna('NoBCR_C')
    else:
        var_df['IGHisotype'] = 'NoBCR'

    if focus_cdr3:
        var_df = var_df.groupby(['CDR3', 'Chain', 'IGHisotype'])['n_cells'].sum().reset_index().rename({'CDR3': 'ID'}, axis=1)
    else:
        var_df['ID'] = var_df.index
        var_df = var_df[['ID', 'CDR3' ,'Chain', 'IGHisotype', 'n_cells']].reset_index(drop=True).copy()
    var_df = var_df[var_df['n_cells'] != 0].copy()
    return var_df


def GroupbyCheck(
    adata,
    groupby,
    groups = None,
    allow_single = True,
    only_ck_groupby = False
):
    if groupby is None:
        raise ValueError("The 'groupby' parameter must be specified and cannot be None.")
    if groupby not in adata.obs.columns:
        valid_list = list(adata.obs.columns.unique())
        raise ValueError(f"'{groupby}' does not exist in adata.obs_names.Only valid in {valid_list}")
    if not only_ck_groupby:
        unique_groups = list(adata.obs[groupby].unique())
        if groups is None:
            candidate = unique_groups
        elif isinstance(groups, list):
            LiteralCheck(groups, unique_groups)
            candidate = groups
        else:
            if not allow_single:
                raise ValueError("The 'groups' parameter must be a list of two or more group names.")
            else:
                LiteralCheck(groups, unique_groups)
                candidate = [groups]
        return candidate