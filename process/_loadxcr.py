from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from typing import Literal, Union, Optional

from .._utils import LiteralCheck



def load_xcr(
    clone_file: str,
    cloneId_barcode_file: str,
    focus_functional_clone: bool = True,
    focus_chain: Optional[Union[Literal['TCR', 'BCR']]] = None,
    from_mixcr: bool = True,
    add_mutation: bool = True,
) -> AnnData:
    """
    Bulid a AnnData object from MIXCR report.

    Args:
        mixcr_file: str
            File path to mixcr contig file. 

        cloneId_barcode_file: str
            File path to clone-barcode file. (split by '\\t')

        focus_functional_clone: Union[Literal['TCR', 'BCR'], None]
            Whether to only focus on functional clones. If True, will filter all clones without functionality.

        from_mixcr: bool
            If clone_file comes from mixcr assemble, it is True. If it is a customized clone_file, The first three 
            columns must be 'cloneId', 'topChains' and 'aaSeqCDR3', and the other columns are optional.

        add_mutation: bool
            If True, mutation information will be extracted from the clone file. If from_mixcr is True, it will be 
            extracted based on the mixcr result. If from_mixcr is False, additional V_align, J_align, and nSeqCDR3 
            need to be provided in the clone file.
        
    Returns:
        AnnData: A AnnData of spatial immune repertoire.
    """
    
    validlist = ['TCR','BCR']
    LiteralCheck(focus_chain, validlist = validlist)

    clone_df = _read_clone_file(
        clone_file,
        focus_functional_clone = focus_functional_clone,
        focus_chain = focus_chain,
        from_mixcr = from_mixcr,
    )
    clone_df = clone_df.reset_index(drop = True)
    cloneId_barcode_df = _read_cloneId_barcode_file(cloneId_barcode_file)

    # ensure var_names is not numeric
    clone_is_numeric = pd.api.types.is_numeric_dtype(cloneId_barcode_df['cloneId'])
    if clone_is_numeric:
        cloneId_barcode_df['cloneId'] = 'clone_' + cloneId_barcode_df['cloneId'].map(str)
        clone_df['cloneId'] = 'clone_' + clone_df['cloneId'].map(str)
    
    # intersection cloneId_barcode_df and clone_df
    cloneId_barcode_df = cloneId_barcode_df[cloneId_barcode_df['cloneId'].isin(clone_df['cloneId'])].copy()
    if cloneId_barcode_df.shape[0] == 0:
        raise ValueError('Not found any overlap between cloneId_barcode and mixcr.')

    # bulid adata
    adata = _bulid_adata(cloneId_barcode_df)

    aaSeqCDR3_dict = dict(zip(clone_df['cloneId'],clone_df['aaSeqCDR3']))
    adata.var['CDR3'] = adata.var_names.map(aaSeqCDR3_dict)
    adata.var['CDR3_lenth'] = adata.var['CDR3'].map(len)

    Func_dict = dict(zip(clone_df['cloneId'],clone_df['Func']))
    adata.var['Func'] = adata.var_names.map(Func_dict)

    topChains_dict = dict(zip(clone_df['cloneId'],clone_df['topChains']))
    adata.var['Chain'] = adata.var_names.map(topChains_dict)
    if focus_chain != 'TCR' and 'isotype' in clone_df.columns and 'IGH' in adata.var['Chain'].unique():
        isotype_dict = dict(zip(clone_df['cloneId'],clone_df['isotype']))
        adata.var['IGHisotype'] = adata.var_names.map(isotype_dict)
        adata.var['IGHisotype'] = adata.var['IGHisotype'].map(
            {
                'A':'IGHA','H':'IGHG','D':'IGHD','E':'IGHE','M':'IGHM',
                'a':'IGHA','h':'IGHG','d':'IGHD','e':'IGHE','m':'IGHM',
                'igA':'IGHA','igH':'IGHG','igD':'IGHD','igE':'IGHE','igM':'IGHM',
                'iga':'IGHA','igh':'IGHG','igd':'IGHD','ige':'IGHE','igm':'IGHM',
                'IgA':'IGHA','IgH':'IGHG','IgD':'IGHD','IgE':'IGHE','IgM':'IGHM',
                'IGA':'IGHA','IGH':'IGHG','IGD':'IGHD','IGE':'IGHE','IGM':'IGHM',
                'IGHA':'IGHA','IGHG':'IGHG','IGHD':'IGHD','IGHE':'IGHE','IGHM':'IGHM',
            }
            )
        # adata.var['IGHisotype'] = adata.var['IGHisotype'].fillna('None')
    adata.uns['immu_repertoire'] = {}
    adata.uns['immu_repertoire']['args'] = {}
    adata.uns['immu_repertoire']['args']['from_mixcr'] = from_mixcr

    # add mutation info:
    adata.uns['immu_repertoire']['args']['mutated'] = False
    if add_mutation:
        if from_mixcr:
            if 'allVHitsWithScore' not in clone_df.columns or 'allJHitsWithScore' not in clone_df.columns:
                raise KeyError("Not Found 'allVHitsWithScore' or 'allJHitsWithScore' in clone file")
            else:
                mutation = pd.DataFrame(
                    clone_df.apply(lambda x : _bulid_mutation_df(x),axis = 1).tolist()
                    )
                clone_df = pd.concat([clone_df,mutation],axis = 1)
                clone_df['mutation_status'] = (clone_df['Mutation'] == 0)
                clone_df['mutation_status'] = clone_df['mutation_status'].map({True:'germline',False:'mutated'})
                adata.uns['immu_repertoire']['args']['mutated'] = True
                
        else:
            columns_needs = ['V_align','J_align','nSeqCDR3','mutation_status']
            for columns_need in columns_needs:
                if columns_need not in clone_df.columns:
                    raise KeyError(f"Not Found {columns_need} in clone file")
            adata.uns['immu_repertoire']['args']['mutated'] = True

            
        
    clone_df = clone_df.fillna('None')
    if from_mixcr:
        clone_df['clonename'] = clone_df['aaSeqCDR3']+'@'+clone_df['V_align']+'@'+clone_df['J_align']
    adata.uns['immu_repertoire']['clone_info'] = clone_df.copy()
    if add_mutation:
        mutation_dict = dict(zip(clone_df['cloneId'],clone_df['mutation_status']))
        adata.var['mutation'] = adata.var_names.map(mutation_dict)
    return adata.copy()


def _getchain(
    clonotypes_df: pd.DataFrame
):
    print("> Not Found \"topChains\" infomation from MIXCR. Get topChains from \"allVHitsWithScore\" ")
    clonotypes_df['topChains'] = clonotypes_df['allVHitsWithScore'].apply(lambda x : x[:3])
    return clonotypes_df

def _getisotype(x):
    if isinstance(x,float):
        return np.nan
    elif x.startswith('IGH'):
        return x[:4]
    else:
        return np.nan
        
def _get_isotype_from_mixcr(
    clonotypes_df: pd.DataFrame
):
    if 'isotypePrimary' not in clonotypes_df:
        clonotypes_df['isotype'] = clonotypes_df['allCHitsWithScore'].map(_getisotype)
    else:
        clonotypes_df['isotype'] = clonotypes_df['isotypePrimary'].map({'IgA':'IGHA','IgE':'IGHE','IgD':'IGHD','IgM':'IGHM','IgG':'IGHG'})
    return clonotypes_df

def _read_clone_file(
    clone_file : str,
    focus_functional_clone: bool = True,
    focus_chain: Union[None,Literal['TCR', 'BCR']] = None,
    from_mixcr: bool = True,
):

    clone_df = pd.read_csv(clone_file,sep = '\t')
    if "allVHitsWithScore" not in clone_df.columns:
        print(f'Warming: {clone_file} may not from Mixcr assemble.')
    # check mixcr format
    if 'Unnamed: 0' in clone_df.columns:
        clone_df = clone_df.set_index('Unnamed: 0')
        clone_df.index.name = None

    if 'topChains' not in clone_df.columns and from_mixcr:
        clone_df = _getchain(clone_df)
    elif 'topChains' not in clone_df.columns and not from_mixcr:
        raise KeyError(f'Not found \'topChains\' in clone_file columns. Only Found {clone_df.columns}')
        
    clone_df['Func'] = clone_df['aaSeqCDR3'].map(lambda x : 'Lose Function' if '*' in x or '_' in x else 'Functional')
    
    if focus_functional_clone:
        clone_df = clone_df[clone_df['Func'] == 'Functional'].copy()
    
    if focus_chain!= None:
        if focus_chain == 'BCR':
            clone_df = clone_df[clone_df['topChains'].isin(['IGH','IGL','IGK'])].copy()
        elif focus_chain == 'TCR':
            clone_df = clone_df[clone_df['topChains'].isin(['TRA','TRB'])].copy()
        else:
            raise ValueError(f"Focus_chain must be BCR or TCR. Invalid parameter {focus_chain}")
    else:
        clone_df = clone_df[clone_df['topChains'].isin(['TRA','TRB','IGH','IGL','IGK'])].copy()
        
    if focus_chain != 'TCR' and from_mixcr:
        clone_df = _get_isotype_from_mixcr(clone_df)
    
    return clone_df


def _read_cloneId_barcode_file(
    cloneId_barcode_file: str
):
    raw_cloneId_barcode_df = pd.read_csv(cloneId_barcode_file,header=None,sep = '\t')
    columns_lenth = raw_cloneId_barcode_df.shape[1]
    # check cloneId_barcode format
    if columns_lenth == 2:
        raw_cloneId_barcode_df.columns = ['barcode','cloneId']
        raw_cloneId_barcode_df['count'] = 1
    elif columns_lenth == 3:
        raw_cloneId_barcode_df.columns = ['barcode','cloneId','count']    
    is_numeric = pd.api.types.is_numeric_dtype(raw_cloneId_barcode_df['count'])
    if not is_numeric:
        raise KeyError("Found Valid dtype in cloneId_barcode_file. ")
    return raw_cloneId_barcode_df

def _bulid_adata(
    inputdf: pd.DataFrame
):
    # inputdf format: barcode cloneId count [sep = '\t']
    col1,col2,col3 = inputdf.columns
    row = inputdf[col1].astype('category').cat.codes
    col = inputdf[col2].astype('category').cat.codes
    data = inputdf[col3].values
    coo = coo_matrix((data, (row, col)))
    csr = coo.tocsr()
    adata = AnnData(X = csr)
    adata.obs_names = inputdf[col1].astype('category').cat.categories
    adata.var_names = inputdf[col2].astype('category').cat.categories

    return adata.copy()


# mutation

def _mutation_extract(
    Alignments_str: str,
):
    queryLength_list = []
    mutation_freq_list = []
    alignmentScore_list = []
    Alignments = Alignments_str.split(';')[0]
    info = Alignments.split('|')
    queryLength = int(info[4]) - int(info[3])
    mutations = info[5]
    alignmentScore = float(info[6])
    subs = mutations.count('S')
    dels =  mutations.count('D')
    ins = mutations.count('I')
    mutation_freq = subs + dels + ins
    queryLength_list.append(queryLength)
    mutation_freq_list.append(mutation_freq)
    alignmentScore_list.append(alignmentScore)
    queryLength_sum = sum(queryLength_list)
    mutation_freq_sum = sum(mutation_freq_list)
    if queryLength_sum != 0:
        mutation_perc = mutation_freq_sum / queryLength_sum
    else:
        mutation_perc = 0
    return queryLength_sum,mutation_freq_sum,mutation_perc,np.median(alignmentScore_list),mutations


def _bulid_mutation_df(
    row: pd.Series
):
    v = row['allVHitsWithScore']
    j = row['allJHitsWithScore']

    j_align = row['allJAlignments']
    jql,jmf,jmp,js,jmutations = _mutation_extract(j_align)

    v_align = row['allVAlignments']
    vql,vmf,vmp,vs,vmutations = _mutation_extract(v_align)

    # not for mutation but for align gene
    d_align = row['allDHitsWithScore']
    if not isinstance(d_align,float):
        d_align  = d_align.split(';')[0].split('(')[0]
    c_align = row['allCHitsWithScore']
    if not isinstance(c_align,float):
        c_align  = c_align.split(';')[0].split('(')[0]
    
    if jmutations == None:
        jmutations = '-'
    if vmutations == None:
        vmutations = '-'
    total_mutation = jmf + vmf
    tmp = {
        'V_align': v.split(';')[0].split('(')[0],
        'V_mutated': vmutations,
        'V_mutated_count':vmf,
        'J_align': j.split(';')[0].split('(')[0],
        'J_mutation': jmutations,
        'J_mutated_count':jmf,
        'Mutation': total_mutation,
        'D_align': d_align,
        'C_align': c_align,
        }
    return tmp