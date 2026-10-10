from typing_extensions import Literal
import pandas as pd
import numpy as np
from anndata import AnnData
from collections import Counter
from sklearn.neighbors import KDTree,NearestNeighbors
from scipy import stats
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

from typing import Optional, Union, List, Dict


from .._utils import LiteralCheck, extract_var, GroupbyCheck
from ._query import query


########################################
#####  aggregate
########################################

def aggregate(adata,groups,groupby,method=None,k = None,d = None,cutoffmode = 'mean',basis = 'spatial',aggregate_name = 'aggregate'):
    if method not in ['knn', 'radius']:
        raise ValueError("method must be either 'knn' for k-nearest neighbors or 'radius'")
    candidate = GroupbyCheck(adata, groupby, groups = groups,allow_single = True)

    expandeadata = adata[adata.obs[groupby].isin(candidate)]
    clonedf = expandeadata.obs.copy()
    if 'x' in clonedf.columns and 'y' in clonedf.columns:
        cdr3loc = clonedf[['x','y']].values
    else:
        if basis in expandeadata.obsm:
            cdr3loc = expandeadata.obsm[basis]
        elif  'X_'+basis in expandeadata.obsm :
            cdr3loc = expandeadata.obsm['X_'+basis]
        else:
            raise ValueError(f"The basis '{basis}' is not found in adata.obsm")

    tree = KDTree(cdr3loc)
    if method == 'knn' and k is not None:
        query = tree.query(cdr3loc, k = k)
        meandistance = query[0].mean(axis = 1)
        if cutoffmode == 'median':
            clonedf['clean'] = (meandistance <= np.median(meandistance)).reshape(-1, 1)
        elif cutoffmode == 'mean':
            clonedf['clean'] = (meandistance <= np.mean(meandistance)).reshape(-1, 1)
        elif cutoffmode.startswith('p'):
            p_num = cutoffmode.split('p')[1]
            if p_num.isdigit():
                p_num = int(p_num)
            else:
                raise KeyError(f"If cutoffmode starts with 'p', it must be followed by an integer.")
            clonedf['clean'] = (meandistance <= np.percentile(meandistance,p_num)).reshape(-1, 1)

    elif method == 'radius' and d is not None:
        nn = NearestNeighbors(radius = d)
        nn.fit(cdr3loc)
        indices = nn.radius_neighbors(cdr3loc, return_distance=False)
        query = np.array([len(index) for index in indices])
        if cutoffmode == 'median':
            mean_spot_count = np.mean(query)
            clonedf['clean'] = (query > mean_spot_count).reshape(-1, 1)
        elif cutoffmode == 'mean':
            median_spot_count = np.median(query)
            clonedf['clean'] = (query > median_spot_count).reshape(-1, 1)
        elif cutoffmode.startswith('p'):
            p_num = cutoffmode.split('p')[1]
            if p_num.isdigit():
                p_num = int(p_num)
            else:
                raise KeyError(f"If cutoffmode starts with 'p', it must be followed by an integer.")
            p_spot_count = np.percentile((query),p_num)
            clonedf['clean'] = (query >= p_spot_count).reshape(-1, 1)
    else:
        raise ValueError("Appropriate 'k' or 'd' value must be provided based on the selected mode.")

    adata.obs[aggregate_name] = (adata.obs_names.isin(clonedf[clonedf['clean']].index))
    adata.obs[aggregate_name] = adata.obs[aggregate_name].map({True:aggregate_name})



########################################
#####  intersection
########################################
def intersection(
    adata: AnnData, 
    groupby: str, 
    chain: Optional[str] = None, 
    focus_cdr3: bool = False, 
    groups: Optional[Union[str, List[str]]] = None, 
    mode: Literal['inner', 'outer', 'total'] = 'outer'
) -> pd.DataFrame:
    """
    Extract CDR3 intersections across specified groups.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    
    groupby : str
        Column name in `adata.obs` used for grouping cells (e.g., clusters, annotations).
    
    chain : Optional[str], optional
        Specifies the chain to focus on. If None, all chains are considered.
    
    focus_cdr3 : bool, optional
        Whether to focus on CDR3 sequences. If False, the clone is considered.
    
    groups : Optional[Union[str, List[str]]], optional
        Groups to consider for comparison. If None, all groups in `groupby` are used.
    
    mode : Literal['inner', 'outer', 'total'], optional
        Mode of intersection. Determines whether to return inner, outer, or total intersections.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing the intersection results based on the chosen mode.
    """
    
    candidate = GroupbyCheck(adata, groupby = groupby, groups = groups)
    
    cdr3_dict = {}
    
    for i in candidate:
        cells = adata[adata.obs[groupby] == i].obs_names
        region_cdr3s = query(adata, obs_names=cells, obs_chain=chain, focus_cdr3=focus_cdr3, flatten=True)
        cdr3_dict[i] = set(region_cdr3s)
    df = _extract_intersection(cdr3_dict, mode=mode)
    
    return df


def _extract_intersection(
    cdr3_dict: Dict,
    mode: Literal['inner', 'outer', 'total'],
):
    if mode == 'outer':
        unique_elements = {}
        element_counts = {}

        for key, elements in cdr3_dict.items():
            for element in elements:
                if element in element_counts:
                    element_counts[element].add(key)
                else:
                    element_counts[element] = {key}
        for element, keys in element_counts.items():
            if len(keys) == 1:
                unique_elements[element] = list(keys)[0]
        result_df = pd.DataFrame(list(unique_elements.items()), columns=['ID', 'Occurance'])
    elif mode == 'inner':
        common_elements = set.intersection(*cdr3_dict.values())
        result_df = pd.DataFrame(list(common_elements), columns=['CDR3'])
        result_df['Occurance'] = 'ALL'
    elif mode == 'total':
        element_categories = {}
        all_categories = set(cdr3_dict.keys())

        for category, elements in cdr3_dict.items():
            for element in elements:
                if element not in element_categories:
                    element_categories[element] = set()
                element_categories[element].add(category)
        for element, categories in element_categories.items():
            if categories == all_categories:
                element_categories[element] = 'ALL'
            else:
                sorted_categories = '+'.join(sorted(categories))
                element_categories[element] = sorted_categories
        result_df = pd.DataFrame(list(element_categories.items()), columns=['ID', 'Occurance'])
    else:
        raise KeyError("Mode must be 'outer' or 'inner' or 'total' ")
    return result_df


########################################
#####  diversity
########################################

from typing import Optional, Literal
import pandas as pd
from anndata import AnnData

def diversity(
    adata: AnnData,
    groupby: str,
    chain: Optional[str] = None,
    focus_cdr3: bool = False,
    method: Literal['gini', 'shannon', 'simpson'] = 'shannon'
) -> pd.DataFrame:
    """
    Compute diversity indices (Shannon, Gini, or Simpson) for groups within an AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with clone information.
    
    groupby : str
        Column name in `adata.obs` to group cells by (e.g., clusters or annotations).
    
    chain : Optional[str], optional
        Specifies the chain to focus on. If None, all chains are considered.
    
    focus_cdr3 : bool, optional
        Whether to focus on CDR3 sequences or clones.
    
    method : Literal['gini', 'shannon', 'simpson'], optional
        The diversity index to compute. Default is 'shannon'.
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing the computed diversity values for each group.
    """
    
    candidate = GroupbyCheck(adata, groupby = groupby)
    
    result_list = []
    clonedict = query(
            adata, 
            obs_names=adata.obs_names, 
            obs_chain=chain, 
            focus_cdr3=focus_cdr3
        )

    for group in candidate:
        subadata = adata[adata.obs[groupby] == group]
        groupcell = subadata.obs_names
        total_clones = []

        for cell in groupcell:
            clones = clonedict.get(cell,0)
            if len(clones) == 0:
                continue
            else:
                total_clones.extend(clones)
        if len(total_clones) != 0:
            diversity_values = _diversity(total_clones, method=method)
        else:
            diversity_values = np.nan
        result_list.append([group, diversity_values])
    
    result = pd.DataFrame(result_list, columns=['group', method]).sort_values(method, ignore_index=True)
    return result


def _diversity(
    labels: Union[pd.DataFrame, pd.Series, List, np.ndarray],
    method: Literal['gini','shannon','simpson']
):
    
    if isinstance(labels, pd.DataFrame):  #  pd.DataFrame means a:30,b:50
        if labels.shape[1] != 2:
            raise ValueError("DataFrame must have exactly two columns: one for labels and one for counts.")
        label_counts = dict(zip(labels.iloc[:, 0], labels.iloc[:, 1]))
    elif isinstance(labels, pd.Series):  # other means  a,a,a,a,a····,b,b,b,b 
        if pd.api.types.is_numeric_dtype(labels.dtype):
            label_counts = labels.to_dict()
        else:
            label_counts = labels.value_counts().to_dict()
    elif isinstance(labels, (list, np.ndarray)):
        label_counts = Counter(labels)
    else:
        raise TypeError("Input must be a list, numpy array, pandas Series, or pandas DataFrame.")
    
    if method.lower() == 'gini':
        return _calc_gini(label_counts)
    elif method.lower() == 'shannon':
        return _calc_shannon(label_counts)
    elif method.lower() == 'simpson':
        return _calc_simpson(label_counts)
    else:
        raise ValueError("Unsupported method. Choose 'gini', 'shannon', or 'simpson'.")


def _calc_gini(label_counts):
    x = list(label_counts.values())
    mad = np.abs(np.subtract.outer(x, x)).mean()
    rmad = mad/np.mean(x)
    g = 0.5 * rmad
    return g

def _calc_shannon(label_counts):
    num = sum(label_counts.values())
    if num == 0:
        return 0
    shannon = -sum((count / num) * np.log(count / num) for count in label_counts.values())
    return shannon

def _calc_simpson(label_counts):
    num = sum(label_counts.values())
    if num == 0:
        return 0
    simpson = sum((count / num) ** 2 for count in label_counts.values())
    return simpson

###########
## smooth
##########

def get_min_distance(adata,basis = 'spatial'):
    locations = adata.obsm[basis].copy()
    distances = np.linalg.norm(locations[:, np.newaxis] - locations, axis=2)
    np.fill_diagonal(distances, np.inf)
    min_distances = np.min(distances, axis=1)
    average_min_distance = np.mean(min_distances)
    return average_min_distance

def dist_decay(data, kernel='Gaussion', alpha=1):
    r = data.mean() / alpha

    if kernel == 'Gaussian':
        data = np.exp(-data ** 2 / (2 * r ** 2))

    if kernel == 'Exponential':
        data = np.exp(-data / (2 * r ** 2))

    if kernel == 'Laplacian':
        data = np.exp(-data / (2 * r))

    if kernel == 'Cauchy':
        data = np.exp(1 / (data ** 2 / r + 1))

    if kernel == 'distance':
        data = 1 / (1 + data)

    if kernel == 'distance2':
        data = 1 / (1 + data ** 2)
    return data

def smooth(
    adata: AnnData,
    obs_name: str,
    alpha: float = 1,
    kernel: Optional[str] = 'Gaussian',
    smooth_distance: Optional[float] = None,
    min_round_spot: int = 1,
    fillna: bool = False,
    basis: str = 'spatial',
    include_self: bool = True,
    add_obs: str = 'smooth',
    denoise: bool = False,
    denoise_k: int = 10,
    denoise_cutoffmode: str = 'mean',
    norm_by: Optional[str] = None,
) -> None:
    """
    Smooths data by calculating a weighted average of neighboring points based on distance decay.

    Parameters:
    ----------
    adata : AnnData
        An AnnData object with spatial coordinates and the target observation values, commonly used for single-cell spatial data.
        
    obs_name : str
        The name of the column in `adata.obs` containing the values to smooth.
        
    alpha : float, optional, default=1
        The distance decay coefficient that controls the rate of decay. Higher values of alpha cause smoothing to focus more on closer neighbors.
        
    kernel : Optional[str], optional, default='Gaussian'
        The type of distance decay function to apply. Options include:
        - 'Gaussian': Gaussian decay
        - 'Exponential': Exponential decay
        - 'Laplacian': Laplacian decay
        - 'Cauchy': Cauchy decay
        - 'distance': Decay as 1/(1 + distance)
        - 'distance2': Decay as 1/(1 + distance^2)

    smooth_distance : Optional[float], optional, default=None
        The smoothing radius distance. If None, it defaults to twice the minimum neighbor distance. If provided, it is set to `max(smooth_distance, min_distance)`.
        
    min_round_spot : int, optional, default=5
        The minimum number of neighboring points required to retain a cell. If a cell has fewer neighbors, it is excluded.
        
    basis : str, optional, default='spatial'
        The name of the coordinate basis to use for smoothing, typically defined in `adata.obsm`.

    add_obs : str, optional, default='smooth'
        The name of the column in `adata.obs` where the smoothed feature values will be added.

    Returns:
    ----------
    None : Updates `adata.obs` by adding the smoothed values in a new column specified by `add_obs`.
    """

    GroupbyCheck(adata,groupby = obs_name,only_ck_groupby=True)

    if smooth_distance == None:
        min_distance = get_min_distance(adata,basis = basis )
        smooth_distance = 3 * min_distance
        print(f'Smooth_distance: {smooth_distance:.2f},min_distance : {min_distance:.2f}')


    if denoise:
        def update_agg(row,obs_name):
            if row[f'{obs_name}_agg'] == f'{obs_name}_agg':
                return row[obs_name]
            else:
                return 0

        adata.obs[f'{obs_name}_validspot'] = adata.obs[obs_name] > np.mean(adata.obs[obs_name])
        aggregate(adata,groups = True,groupby = f'{obs_name}_validspot',method = 'knn',k = denoise_k,cutoffmode = denoise_cutoffmode,aggregate_name=f'{obs_name}_agg')
        min_value = adata.obs[obs_name].min()
        max_value = adata.obs[obs_name].max()
        adata.obs[f'denoise_{obs_name}'] = (adata.obs[obs_name] - min_value) / (max_value - min_value)
        adata.obs[f'denoise_{obs_name}'] = adata.obs.apply(lambda x: update_agg(x,obs_name), axis = 1)
        del adata.obs[f'{obs_name}_validspot']
        del adata.obs[f'{obs_name}_agg']
        obs_name = f'denoise_{obs_name}'


    if norm_by is None:
        obsdf = adata.obs[obs_name].to_frame()
    else:
        GroupbyCheck(adata,groupby = norm_by,only_ck_groupby = True)
        obsdf = adata.obs[[obs_name,norm_by]].copy()
    obsdf['x'] = adata.obsm[basis][:,0]
    obsdf['y'] = adata.obsm[basis][:,1]
    if fillna:
        obsdf[obs_name] = obsdf[obs_name].fillna(0)
    else:
        obsdf = obsdf.dropna()

    location = obsdf[['x','y']]
    obsdf['rank'] = range(obsdf.shape[0])
    obsdict = dict(zip(obsdf['rank'],obsdf[obs_name]))
    if norm_by is not None:
        obsdf[norm_by] = obsdf[norm_by].fillna(0)
        normdict = dict(zip(obsdf['rank'],obsdf[norm_by]))
    tree = KDTree(location)
    query = tree.query_radius(location, r = smooth_distance,return_distance = True)
    idxs = query[0]
    distances = query[1]

    round_counts = []
    if not include_self:
        distances_dropzero = []
        for distance in distances:
            distance = distance[distance != 0]
            distances_dropzero.append(distance)
            round_counts.append(len(distance))

        idxs_dropzero = []
        for i,idx in enumerate(idxs):
            idx = idx[idx != i]
            idxs_dropzero.append(idx)
    else:
        for distance in distances:
            round_counts.append(len(distance))
        distances_dropzero = distances
        idxs_dropzero = idxs
        
    smooth_list = []
    for now_idx, (distance, idx) in enumerate(zip(distances_dropzero, idxs_dropzero)):
        if norm_by is None:
            if len(idx) != 0:
                weights = dist_decay(distance, kernel = kernel, alpha=1)
                mean_values = np.average([obsdict[i] for i in idx], weights=weights)
            else:
                mean_values = 0  
            smooth_list.append(mean_values)
        else:
            if len(idx) != 0:
                total_values = np.sum([obsdict[i] for i in idx])
                total_background_values = np.sum([normdict[i] for i in idx])
                if total_background_values != 0:
                    norm_values = total_values/total_background_values
                else:
                    norm_values = 0
            else:
                norm_values = 0
            smooth_list.append(norm_values)

    rounddf = pd.DataFrame({'cell':location.index,'roundcell':round_counts})
    obsdf['smooth'] = smooth_list
    obsdf = obsdf.loc[rounddf[rounddf['roundcell'] >= min_round_spot]['cell']]
    resultdict = dict(zip(obsdf.index,obsdf['smooth']))
    adata.obs[add_obs] = adata.obs_names.map(resultdict)
    if fillna:
        adata.obs[add_obs] = adata.obs[add_obs].fillna(0)


########################################
#####  after smooth, rank 
########################################
from sklearn.cluster import KMeans
def kmean_cluster(
    adata: AnnData,
    obs: str,
    cluster:int = 2):
    data = adata.obs[obs].dropna().sort_values()
    arrdata = data.values
    arrdata = arrdata.reshape(-1,1)
    kmeans = KMeans(n_clusters=cluster, random_state = 61).fit(arrdata)
    labels = kmeans.labels_
    cluster = pd.DataFrame({'spot':data.index,'value':data.values,'label':labels})
    cluster['rank'] = range(cluster.shape[0])
    cluster['label'] = cluster['label'].map(str)
    return cluster

def obsrank(
    adata: AnnData,
    obs_name: str,
    n_cluster: int = 2,
    select_top: bool = True,
    add_obs: Optional['str'] = None):
    
    GroupbyCheck(adata,groupby = obs_name,only_ck_groupby=True)

    if add_obs is None:
        add_obs = f"{obs_name}_rank"
    clusterdf = kmean_cluster(adata = adata,obs = obs_name,cluster = n_cluster)
    rank = clusterdf.groupby('label')['value'].mean().sort_values(ascending=False).reset_index()
    rankdict = dict(zip(rank['label'],rank.index))
    clusterdf['label'] = clusterdf['label'].map(rankdict)
    clusterdf['label'] = clusterdf['label'].astype(str)
    if select_top:
        adata.obs[add_obs] = adata.obs_names.isin(clusterdf[clusterdf['label'].isin(['0'])]['spot'])
    else:
        adata.obs[add_obs] = clusterdf.set_index('spot')['label']
    if 'aggragetion' not in adata.uns.keys():
        arg_dict = {}
        adata.uns['aggragetion'] = arg_dict

    arg_dict = {}
    arg_dict.update(adata.uns['aggragetion'])
    arg_dict[add_obs] = {}
    arg_dict[add_obs]['value_from'] = obs_name
    arg_dict[add_obs]['kmean'] = clusterdf
    adata.uns['aggragetion'] = arg_dict
    print(f"Add .. **{add_obs}** in adata.uns['aggragetion'] from obs: {obs_name})")

############
#### density
###########
from collections import defaultdict

def density(
    adata,
    clones = None,
    cdr3s = None, 
    chain = None,
    consider_umi = False,
    density_name = "density", 
    smooth_density = False, 
    smooth_alpha = 1,
    smooth_distance = None,
    norm_by = None,
    result_norm = None,
):
    
    if clones is not None:
        cellwithclones = query(adata,var_names = clones,flatten = False,reverse=True)
    elif cdr3s is not None:
        cellwithclones = query(adata,cdr3_names = cdr3s,flatten = False,reverse=True)
    elif chain is not None:
        vardf = extract_var(adata, chain = chain, focus_cdr3 = False)
        cellwithclones = query(adata,var_names = vardf['ID'],reverse=True)


    cellwithclonescounts = {}
    if not consider_umi:
        for cell,clones in cellwithclones.items():
            cellwithclonescounts[cell] = len(clones)
    else:
        if cdr3s is None:
            for cell,clones in cellwithclones.items():
                cell_umi_count = 0
                for clone in clones:
                    cell_umi_count += adata[cell,clone].X.data[0]
                cellwithclonescounts[cell] = cell_umi_count
        else:
            vardf = extract_var(adata, chain = chain, focus_cdr3 = False)
            cdr32clone = vardf[['CDR3','ID']].groupby('CDR3')['ID'].apply(list).to_dict()
            for cell,cdr3s in cellwithclones.items():
                cell_umi_count = 0
                for cdr3 in cdr3s:
                    for clone in cdr32clone[cdr3]:
                        umi_count = adata[cell,clone].X.data  
                        if umi_count.shape[0] !=0:
                            cell_umi_count += umi_count[0]
                cellwithclonescounts[cell] = cell_umi_count
    adata.obs[density_name] = adata.obs_names.map(cellwithclonescounts)
    if smooth_density:
        smooth(adata,obs_name = density_name, alpha = smooth_alpha, smooth_distance = smooth_distance,fillna = True,min_round_spot = 0,add_obs = density_name,include_self=True, norm_by = norm_by)
    if result_norm:
        df = adata.obs.copy()
        df['col_normalized'] = (df[density_name] - df[density_name].min()) / (df[density_name].max() - df[density_name].min())
        adata.obs[density_name] = df['col_normalized']
        
        
        
###########################
####### pair clone
###########################


def _perm_worker(args):
    pool_a, pool_b, n_a, n_b, r, n_perm, seed, name_a, name_b, obs_count = args
    rng = np.random.default_rng(seed)
    perm_counts = np.zeros(n_perm, dtype=int)

    for i in range(n_perm):
        idx_a = rng.choice(len(pool_a), size=n_a, replace=False)
        idx_b = rng.choice(len(pool_b), size=n_b, replace=False)
        sample_a = pool_a[idx_a]
        sample_b = pool_b[idx_b]

        tree = KDTree(sample_b)
        neighbors = tree.query_radius(sample_a, r=r)
        perm_counts[i] = sum(1 for n in neighbors if len(n) > 0)

    mu = perm_counts.mean()
    std = perm_counts.std(ddof=1) if n_perm > 1 else 0.0

    if std > 1e-10:
        z = (obs_count - mu) / std
        p_value = 2 * (1 - stats.norm.cdf(abs(z)))
        p_value = min(p_value, 1.0)
    else:
        p_value = 0.0 if obs_count != mu else 1.0
    z_score = (obs_count - mu) / std if std > 1e-10 else np.inf

    return {
        'name_a': name_a,
        'name_b': name_b,
        'observed': obs_count,
        'expected': mu,
        'std': std,
        'p_value': p_value,
        'z_score': z_score,
        'n_a': n_a,
        'n_b': n_b
    }


def spatial_pairing(
    df_a,
    df_b,
    r=5,
    n_perm=1000,
    random_state=None,
    n_jobs=-1,
    chunksize=None,
    progress_desc='Permutation test',
    stat=True,
    min_observed=0
):

    df_a = df_a[['name', 'x', 'y']].drop_duplicates().copy()
    df_b = df_b[['name', 'x', 'y']].drop_duplicates().copy()

    for df in [df_a, df_b]:
        df['x'] = df['x'].astype(float)
        df['y'] = df['y'].astype(float)

    all_coords = pd.concat([df_a[['x', 'y']], df_b[['x', 'y']]]).drop_duplicates().reset_index(drop=True)
    coord_arr = all_coords[['x', 'y']].values
    coord_to_idx = {(coord_arr[i, 0], coord_arr[i, 1]): i for i in range(len(coord_arr))}

    idx_to_as = {}
    for name, x, y in df_a[['name', 'x', 'y']].values:
        idx = coord_to_idx[(x, y)]
        idx_to_as.setdefault(idx, []).append(name)

    idx_to_bs = {}
    for name, x, y in df_b[['name', 'x', 'y']].values:
        idx = coord_to_idx[(x, y)]
        idx_to_bs.setdefault(idx, []).append(name)

    tree = KDTree(coord_arr)
    neighbors_list = tree.query_radius(coord_arr, r=r)

    # ---------- 3. 观察值统计 ----------
    observed = {}
    for idx, a_names in idx_to_as.items():
        for n_idx in neighbors_list[idx]:
            if n_idx not in idx_to_bs:
                continue
            b_names = idx_to_bs[n_idx]
            for a_name in a_names:
                for b_name in b_names:
                    key = (a_name, b_name)
                    observed[key] = observed.get(key, 0) + 1

    # 直接过滤 observed < min_observed 的 pair
    if min_observed > 0:
        observed = {k: v for k, v in observed.items() if v >= min_observed}

    if not observed:
        cols = ['name_a', 'name_b', 'observed', 'expected', 'std', 'p_value', 'z_score', 'n_a', 'n_b']
        return pd.DataFrame(columns=cols)

    n_a_dict = df_a.groupby('name').size().to_dict()
    n_b_dict = df_b.groupby('name').size().to_dict()

    # ---------- 4. 快速模式：跳过置换检验 ----------
    if not stat:
        results = []
        for (name_a, name_b), obs_count in observed.items():
            results.append({
                'name_a': name_a,
                'name_b': name_b,
                'observed': obs_count,
                'expected': np.nan,
                'std': np.nan,
                'p_value': np.nan,
                'z_score': np.nan,
                'n_a': n_a_dict[name_a],
                'n_b': n_b_dict[name_b]
            })
        return pd.DataFrame(results)

    # ---------- 5. 统计模式：多进程置换检验 ----------
    pool_a = df_a[['x', 'y']].drop_duplicates().values.astype(float)
    pool_b = df_b[['x', 'y']].drop_duplicates().values.astype(float)

    rng = np.random.default_rng(random_state)
    tasks = []
    for (name_a, name_b), obs_count in observed.items():
        n_a = n_a_dict[name_a]
        n_b = n_b_dict[name_b]

        if n_a > len(pool_a):
            raise ValueError(f"{name_a} 需要 {n_a} 个坐标，但 A 组背景池仅 {len(pool_a)} 个。")
        if n_b > len(pool_b):
            raise ValueError(f"{name_b} 需要 {n_b} 个坐标，但 B 组背景池仅 {len(pool_b)} 个。")

        seed = int(rng.integers(0, 2**32, dtype=np.uint32))
        tasks.append((
            pool_a, pool_b, n_a, n_b, r, n_perm, seed,
            name_a, name_b, obs_count
        ))

    n_jobs = cpu_count() if n_jobs == -1 else int(n_jobs)
    if chunksize is None:
        chunksize = max(1, len(tasks) // (n_jobs * 4))

    results = []
    with Pool(processes=n_jobs) as pool:
        for res in tqdm(
            pool.imap_unordered(_perm_worker, tasks, chunksize=chunksize),
            total=len(tasks),
            desc=progress_desc,
            unit='pair'
        ):
            results.append(res)

    result = pd.DataFrame(results)
    _, fdr_corrected_p_values, _, _ = multipletests(result['p_value'], alpha=0.05, method='fdr_bh')
    result['adjust_p'] = fdr_corrected_p_values
    result = result.sort_values(['adjust_p', 'observed'], ascending=[True, False])
    return result


def XCR_pairing(
    adata,
    alignment,
    receptor_type='BCR',
    r=5,
    n_perm=1000,
    random_state=42,
    n_jobs=-1,
    chunksize=None,
    stat=True,
    min_observed=2,
    progress_desc="pairing"
):
    import os
    import pandas as pd
    
    if receptor_type == 'BCR':
        chain_map = {'IGK': 'L', 'IGL': 'L', 'IGH': 'H'}
        chain_a, chain_b = 'H', 'L'
    elif receptor_type == 'TCR':
        chain_map = {'TRA': 'L', 'TRB': 'H'}
        chain_a, chain_b = 'H', 'L'
    else:
        raise ValueError("receptor_type must be 'BCR' or 'TCR'")
    
    if 'immu_repertoire' not in adata.uns or 'clone_info' not in adata.uns['immu_repertoire']:
        raise KeyError("adata.uns['immu_repertoire']['clone_info'] not found")
    
    clonedf = adata.uns['immu_repertoire']['clone_info']
    
    df = alignment.copy()
    required = {'x', 'y', 'cloneId'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"align_file missing required columns: {missing}")
    
    clonechains = dict(zip(clonedf['cloneId'], clonedf['topChains'].map(chain_map)))
    df['cloneId'] = 'clone_' + df['cloneId'].astype(str)
    df = df[df['cloneId'].isin(clonedf['cloneId'])].copy()
    df['x'] = df['x'].astype(int)
    df['y'] = df['y'].astype(int)
    df['chain'] = df['cloneId'].map(clonechains)
    
    df_a = df[df['chain'] == chain_a][['cloneId', 'x', 'y']].rename(columns={'cloneId': 'name'})
    df_b = df[df['chain'] == chain_b][['cloneId', 'x', 'y']].rename(columns={'cloneId': 'name'})
    
    if df_a.empty or df_b.empty:
        raise ValueError(f"{chain_a} chain or {chain_b} chain data is empty, check receptor_type and align_file")
    
    return spatial_pairing(
        df_a=df_a,
        df_b=df_b,
        r=r,
        n_perm=n_perm,
        random_state=random_state,
        n_jobs=n_jobs,
        chunksize=chunksize,
        progress_desc=progress_desc,
        stat=stat,
        min_observed=min_observed
    )