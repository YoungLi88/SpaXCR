import Levenshtein
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

import pandas as pd
import numpy as np


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import to_hex
from tqdm import tqdm
from multiprocessing import Pool


from typing import Optional, Union, List, Dict, Tuple

from ._utils import _get_color_list,get_colors
from ..extra.radialtree import radialTreee


def CDR3cluster(
    cdr3s: Union[List, pd.DataFrame],
    threshold: Optional[int] = None,
    linkage_method: str = 'ward',
    return_distance: bool = False,
    insert_cost=1,
    replace_cost=1,
    threads = 1,
):
    # input format check
    if isinstance(cdr3s,list):
        cdr3_column_name = 'CDR3'
        cdr3s = pd.DataFrame(cdr3s,columns = ['CDR3'])
    elif isinstance(cdr3s, pd.Series):
        cdr3_column_name = cdr3s.name if cdr3s.name else 'CDR3'
        cdr3s = cdr3s.to_frame(name=cdr3_column_name)
        if cdr3_column_name.lower() != 'cdr3':
            cdr3s.columns = ['CDR3']
            cdr3_column_name = 'CDR3'
        
    elif isinstance(cdr3s,pd.DataFrame):
        if 'CDR3' in cdr3s:
            cdr3_column_name = 'CDR3'
        elif 'cdr3' in cdr3s:
            cdr3_column_name = 'cdr3'
        else:
            raise KeyError(f"If provide pd.DataFrame, DataFrame must contain CDR3 or cdr3 in columns")

    # cluster 
    cdr3s = cdr3s.sort_values(cdr3_column_name)
    d = _distance(cdr3s[cdr3_column_name],insert_cost = insert_cost,replace_cost = replace_cost,threads = threads)
    Z = linkage(squareform(d), linkage_method)
    if threshold != None:
        clusters = fcluster(Z, threshold, criterion='distance')
        cdr3s['cluster'] = clusters
    if return_distance:
        return Z,d,cdr3s
    else:
        return Z,cdr3s


def levenshtein(a, b, insert_cost=1, replace_cost=1):
    dp = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    
    for i in range(len(a) + 1):
        dp[i][0] = i * insert_cost
    for j in range(len(b) + 1):
        dp[0][j] = j * insert_cost

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] 
            else:
                dp[i][j] = min(
                    dp[i - 1][j] + insert_cost,   
                    dp[i][j - 1] + insert_cost,  
                    dp[i - 1][j - 1] + replace_cost 
                )
    return dp[len(a)][len(b)]

# def _distance(cdr3s,insert_cost=1,replace_cost=1):
#     if not isinstance(cdr3s,list):
#         cdr3s = list(cdr3s)
#     n = len(cdr3s)
#     dist_matrix = np.zeros((n, n))
#     if n > 2000:
#         print(f'Found total {n} CDR3s, It may takes long time')
#     for i in tqdm(range(n),total = n):
#         for j in range(i + 1, n):
#             # dist = Levenshtein.distance(cdr3s[i], cdr3s[j])
#             dist = levenshtein(cdr3s[i], cdr3s[j],insert_cost=insert_cost, replace_cost=replace_cost)
#             dist_matrix[i, j] = dist
#             dist_matrix[j, i] = dist
#     return dist_matrix

def calculate_distance(args):
    i, j, cdr3s, insert_cost, replace_cost = args
    dist = levenshtein(cdr3s[i], cdr3s[j], insert_cost=insert_cost, replace_cost=replace_cost)
    return (i, j, dist)

def _distance(cdr3s, insert_cost=1, replace_cost=1, threads=1):
    if not isinstance(cdr3s, list):
        cdr3s = list(cdr3s)
    n = len(cdr3s)
    dist_matrix = np.zeros((n, n))
    
    if n > 2000 and threads == 1 :
        print(f'Found total {n} CDR3s, It may take a long time.Please use multiple threads')

    if threads == 1:
        for i in tqdm(range(n), total=n):
            for j in range(i + 1, n):
                dist = levenshtein(cdr3s[i], cdr3s[j], insert_cost=insert_cost, replace_cost=replace_cost)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist
    else:
        print(f'Begin with {threads} threads')
        tasks = [(i, j, cdr3s, insert_cost, replace_cost) for i in range(n) for j in range(i + 1, n)]
        with Pool(threads) as pool:
            results = list(tqdm(pool.imap(calculate_distance, tasks,chunksize= 1000), total=len(tasks)))
        for i, j, dist in results:
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix

def _similarity(distance_matrix, method = 'gaussian', beta = None,):
    if method not in ['gaussian', 'linear', 'exponential']:
        raise ValueError("Method must be 'gaussian', 'linear', or 'exponential'.")
    if beta is None and method in ['gaussian', 'exponential']:
        beta = 1.0 / np.mean(distance_matrix)

    if method == 'gaussian':
        similarity_matrix = np.exp(-beta * distance_matrix**2)
    elif method == 'linear':
        similarity_matrix = 1 / (1 + distance_matrix)
    elif method == 'exponential':
        similarity_matrix = np.exp(-beta * distance_matrix)
    return similarity_matrix


def similarity_heatmap(
    cdr3s: pd.DataFrame,
    label: Optional[str] = None,
    cmap: str = 'Spectral_r',
    similarity_method: str = 'gaussian',
    beta: Optional[float] = None,
    linkage_method: str = 'ward',
    figsize: Tuple[int, int] = (4, 4),
    label_colors: Optional[Union[str, List[str], Dict[str, str]]] = None,
    label_show: bool = False,
    threshold: Optional[int] = 50,
    return_df = False,
    insert_cost:int = 1,
    replace_cost: int = 1,
    triu: bool = False,
    return_g: bool = False,
    threads: int = 1,
    **kwarg
):
    """
    Generate a similarity heatmap using the specified parameters and CDR3 clustering.

    Parameters
    ----------
    cdr3s : pd.DataFrame
        DataFrame containing CDR3 sequence data.
    
    label : Optional[str], optional
        Column name for the labels to be used for coloring rows in the heatmap. If cluster, it will cluster
    
    cmap : str, optional
        Colormap to use for the heatmap. Default is 'Spectral_r'.
    
    similarity_method : str, optional
        Method for calculating similarity. Can be 'gaussian', 'linear', or 'exponential'.
    
    beta : Optional[float], optional
        Optional beta value for similarity calculation.
    
    linkage_method : str, optional
        Linkage method to use for clustering. Default is 'ward'.
    
    figsize : tuple, optional
        Figure size for the heatmap. Default is (4, 4).
    
    label_colors : Optional[Union[str, List[str], Dict[str, str]]], optional
        Colors for the labels. If None, default coloring will be applied.
    
    label_show: bool, optional
        If label_show, show label legend.
    
    threshold : Optional[int], optional
        Threshold for clustering CDR3 sequences. If None, no threshold is applied.
        
    return_df: bool, optional
        If return_df, return df.
    
    **kwarg : Additional keyword arguments for `sns.clustermap`.
    
    """
    
    Z, d, cdr3s = CDR3cluster(cdr3s, threshold=threshold, linkage_method=linkage_method, return_distance=True,insert_cost = insert_cost,replace_cost=replace_cost,threads = threads)
    s = _similarity(d, method=similarity_method, beta=beta)
    

    if label != None and label not in cdr3s.columns:
        raise KeyError(f'{label} not in DataFrame.columns')
    
    if label != None:
        label_data = cdr3s[label]
        if label != 'cluster':
            unqiue_label = label_data.unique()
            label_color = get_colors(unqiue_label, colors = label_colors)
            label_color = dict(zip(unqiue_label,label_color))
        else:
            if isinstance(label_colors,str):
                cmap_colors = plt.get_cmap(label_colors).colors
                max_cluster = cdr3s['cluster'].max()
                tmpcmap = [to_hex(cmap_colors[(i) % len(cmap_colors)]) for i in range(1,max_cluster+1)]
                label_color = dict(zip(range(1,max_cluster+1),tmpcmap))   #  {1:'xxx',2:'xxxx'}
            elif isinstance(label_colors,dict):
                label_color = label_colors
            label_color = valid_pallete(label_color,cdr3s,pallete = None)
        row_colors = [label_color[i] for i in label_data]
        kwarg['row_colors'] = row_colors
    
    if 'vmax' not in kwarg:
        kwarg['vmax'] = 1
    if 'vmin' not in kwarg:
        kwarg['vmin'] = 0

    g = sns.clustermap(
        s,
        row_linkage=Z,
        col_linkage=Z,
        cmap=cmap,
        figsize=figsize,
        cbar_pos=(1.05, 0.05, 0.02, 0.35),
        xticklabels=False,
        yticklabels=False,
        **kwarg
    )

    if triu:
        mask = np.triu(np.ones_like(s), k = 1)
        values = g.ax_heatmap.collections[0].get_array().reshape(s.shape)
        new_values = np.ma.array(values, mask=mask)
        g.ax_heatmap.collections[0].set_array(new_values)
        g.ax_col_dendrogram.set_visible(False)
        for i in range(s.shape[0]):
            for j in range(s.shape[1]):
                if mask[i, j]: 
                    g.ax_heatmap.add_patch(plt.Rectangle((j, i), 1, 1, 
                                                        facecolor="white", edgecolor="white"))


    if label!= None: 
        sorted_items = sorted(label_color.items()) 
        legend_handles = [mpatches.Patch(color=color, label=tlabel) for tlabel, color in sorted_items]
        if label_show:
            total_color = len(label_color)
            ncol = total_color//3
            g.ax_heatmap.legend(handles=legend_handles,ncol = ncol, bbox_to_anchor=(1.05, 1.05),loc='upper left',title='Cluster', frameon=False)
    
    if return_df:
        return cdr3s
    if return_g:
        return g
    else:
        plt.show()


import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

def similarity_radialtree(
    cdr3s, 
    threshold: Optional[int] = None,
    linkage_method: str = 'ward',
    pallete: Optional[List[str]] = None,
    cluster_in_meta: bool = False,
    metas: Optional[Union[str, List[str]]] = None,
    color_dicts: Optional[Union[None, Dict[str, str], List[Dict[str, str]]]] = None,
    linecolor: str = 'black',
    linewidth: int = 1,
    figsize: Optional[Tuple[int, int]] = None,
    ax: Optional[plt.Axes] = None,
    save: Optional[str] = None,
    alpha: float = 0.7,
    return_df: bool = False,
    no_plot: bool = False,
    insert_cost: int = 1,
    replace_cost: int = 1,
    rangefrom: tuple = (1,10),
    threads: int = 1,
) -> Optional[plt.Axes]:
    """
    Plots a radial tree with colored clusters based on CDR3 data.

    Parameters
    ----------
    cdr3s : pd.DataFrame
        DataFrame containing CDR3 sequence data.
    
    threshold : int, optional
        Threshold for clustering. Default is 50.
    
    linkage_method : str, optional
        Linkage method for hierarchical clustering. Default is 'ward'.
    
    pallete : Optional[List[str]], optional
        List of colors for the clusters. If None, a default palette is used.
    
    metas : Optional[Union[str, List[str]]], optional
        Metadata for labeling the nodes.
    
    color_dicts : Optional[Union[None, Dict[str, str], List[Dict[str, str]]]], optional
        Dictionary or list of dictionaries for color mapping. If None, default colors are used.
    
    linecolor : str, optional
        Color for the lines in the radial tree. Default is 'black'.
    
    linewidth : int, optional
        Line width for the radial tree. Default is 1.
    
    figsize : Optional[Tuple[int, int]], optional
        Figure size. If None, default is (10, 10).
    
    ax : Optional[plt.Axes], optional
        Axes on which to plot the radial tree. If None, a new figure is created.
    
    save : Optional[str], optional
        Path to save the plot. If None, the plot is not saved.
    
    alpha : float, optional
        Alpha transparency for the plot. Default is 0.7.

    return_df: bool, optional
        If return_df, return df.
    
    Returns
    -------
    Optional[plt.Axes]
        Returns the axes with the plot, or None if a new figure is created.
    """
    if threshold is None:
        p = []
        cluster_nums = []
        Z,clusterdf = CDR3cluster(cdr3s,threshold = None,linkage_method=linkage_method, return_distance=False,insert_cost = insert_cost,replace_cost=replace_cost,threads = threads)
        for thrshold in range(rangefrom[0],rangefrom[1]):
            clusters = fcluster(Z, thrshold ,criterion = 'distance')
            clusterdf['cluster'] = clusters
            cluster_num = clusterdf['cluster'].value_counts().shape[0]
            p.append(clusterdf.shape[0]/cluster_num)
            cluster_nums.append(cluster_num)
        tmpdf = pd.DataFrame({'cluster_num':cluster_nums,'p':p}) 
        tmpdf['x'] = (tmpdf['cluster_num'] - tmpdf['p']).map(abs)
        threshold = tmpdf.sort_values('x',ascending=True).index[0]
        threshold += 1
        print(f'Auto set threshold to {threshold} (range from {rangefrom[0]} to {rangefrom[1]}),if you wanna change the range.please use rangefrom = (1,20)')

        clusters = fcluster(Z, threshold ,criterion = 'distance')
        clusterdf['cluster'] = clusters
        cdr3s = clusterdf
    else:
        Z, cdr3s = CDR3cluster(cdr3s, threshold=threshold, linkage_method=linkage_method, return_distance=False,insert_cost = insert_cost,replace_cost=replace_cost,threads = threads)
    if no_plot:
        return cdr3s
    if pallete is None:
        pallete_color = _get_color_list(k=len(cdr3s['cluster'].unique()))
        cgroup = [f'C'+str(x+1) for x in range(len(pallete_color))]
        pallete = dict(zip(cgroup,pallete_color))  # {'C1':'xxx','C2':'xxxx'}
        cluster_pallete = dict(zip([x+1 for x in range(len(pallete))],pallete_color))#  {1:'xxx',2:'xxxx'}
    elif isinstance(pallete,str):
        cmap_colors = plt.get_cmap(pallete).colors
        max_cluster = cdr3s['cluster'].max()
        cmap = [to_hex(cmap_colors[(i) % len(cmap_colors)]) for i in range(1,max_cluster+1)]
        cluster_pallete = dict(zip(range(1,max_cluster+1),cmap))   #  {1:'xxx',2:'xxxx'}
        pallete = dict(zip([ f'C{i}' for i in cluster_pallete.keys()],cluster_pallete.values()))
    else:
        cluster_pallete = pallete  # cluster_pallete    #  {1:'xxx',2:'xxxx'}
        pallete = dict(zip([ f'C{i}' for i in pallete.keys()],pallete.values()))

    pallete,cluster_pallete = valid_pallete(cluster_pallete,cdr3s,pallete = pallete)
    # remove isolate 

    rZ = dendrogram(Z, no_plot=True, color_threshold=threshold)
    

    # Initialize color step tracking
    before = ""
    result_list = []
    step_C = 0

    for now in rZ['color_list']:
        if now == 'C0':
            result_list.append(now)
        else:
            if now != before:
                step_C += 1
            result_list.append(f'C{step_C}')
            before = now
    rZ['color_list'] = result_list
#     return pallete,cluster_pallete,cdr3s,rZ


    plot_colors_dict = {}
    plot_colors_legends = {}

    if cluster_in_meta:
        if metas == None:
            metas = 'cluster'
        else:
            if isinstance(metas, str):
                metas = [metas,'cluster']
            elif isinstance(metas, list):
                metas.append('cluster')

    if metas is not None:
        if isinstance(metas, str):
            metas = [metas]
        if not isinstance(color_dicts, list):
            color_dicts = [color_dicts]*len(metas)
        elif isinstance(color_dicts,list):
            if cluster_in_meta:
                color_dicts.append('tmp')

        for meta, colors_dict in zip(metas, color_dicts):
            if meta == 'cluster':
                colors_dict = cluster_pallete
            
            if not isinstance(colors_dict, dict):
                colors = get_colors(cdr3s[meta].unique(), None)
                colors_dict = dict(zip(cdr3s[meta].unique(), colors))

            plot_colors_dict.update(
                {
                    meta: np.array(cdr3s[meta].map(colors_dict).tolist())
                }
            )
            plot_colors_legends.update(
                {
                    meta: {'colors': np.array(list(colors_dict.values())), 'labels': list(colors_dict.keys())}
                }
            )

    returnax = True
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (10, 10)
        fig, ax = plt.subplots(figsize=figsize)

    # Plot the radial tree
    radialTreee(
        rZ, 
        ax=ax, 
        colorlabels=plot_colors_dict, 
        colorlabels_legend=plot_colors_legends, 
        pallete=pallete, 
        addlabels=False, 
        alpha=alpha, 
        linecolor=linecolor,
        linewidth=linewidth
    )
    
    plt.tight_layout()

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')

    if returnax:
        return ax
    else:
        plt.show()

    if return_df:
        return cdr3s



def valid_pallete(cluster_pallete,cdr3s,pallete = None):
    cluster_valid = list(cluster_pallete.keys())
    invalid_clusters = cdr3s['cluster'].value_counts()[cdr3s['cluster'].value_counts() == 1].index
    if len(invalid_clusters) == 0 :
        new_cluster_pallete = cluster_pallete
    else:
        for invalid_cluster in invalid_clusters:
            cluster_valid.remove(invalid_cluster)
            new_cluster_pallete = dict(zip(cluster_valid,cluster_pallete.values()))
    
    if pallete!= None:
        new_pallete = dict(zip([f'C{i}' for i in range(1,len(new_cluster_pallete)+1)],new_cluster_pallete.values()))

    for invalid_cluster in cdr3s['cluster'].value_counts()[cdr3s['cluster'].value_counts() == 1].index:
        new_cluster_pallete[invalid_cluster] = 'grey'

    if pallete!= None:
        return new_pallete,new_cluster_pallete
    else:
        return new_cluster_pallete