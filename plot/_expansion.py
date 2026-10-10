import matplotlib.pyplot as plt
from anndata import AnnData
import pandas as pd
import numpy as np
from sklearn.neighbors import KDTree
import scanpy as sc 
import seaborn as sns

from typing import Union, List, Optional, Tuple, Dict
from typing_extensions import Literal


from ..tools import query
from ..tools import expansion as expansion_tl
from .._utils import LiteralCheck, extract_var

from ._utils import get_colors

def expansion(
    adata: AnnData, 
    chain: str, 
    focus_cdr3: bool = False, 
    kind: Literal['bar', 'pie'] = 'pie',
    ntop: int = 5,
    label_loc: Optional[Literal['center', 'egde']] = None,
    groups: Union[List[int], int] = (1, 10), 
    labels: Optional[List[str]] = ['Notexpansion', 'Expansion', 'Hyperexpansion'], 
    ax: Optional[plt.Axes] = None, 
    figsize: Optional[Tuple[float, float]] = None, 
    colors: Union[str, List[str], Dict[str, str]] = 'Set2', 
    fontsize: int = 12, 
    save: Optional[str] = None, 
    **kwarg
):
    """
    Plot expansion data for a given AnnData object and chain type, either as a bar chart or pie chart.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing clone information.
        
    chain : str, optional, default='TRB'
        Chain type to focus on.
        
    focus_cdr3 : bool, optional, default=False
        Whether to focus on CDR3 for expansion analysis.
        
    kind : Literal['bar', 'pie'], optional, default='pie'
        Type of plot ('bar' or 'pie').
        
    ntop : int, optional, default=5
        Number of top items to display in the plot. [for bar plot]
    
    label_loc: Literal['center', 'egde'], optional, default = None
        Only for bar plot. Add label_loc on bar plot.

    groups : Union[List[int], int], optional, default=(1, 10)
        Defines the intervals of cell counts to categorize CDR3 or chains. [for pie plot]
        
    labels : List[str], optional, default=['Notexpansion', 'Expansion', 'Hyperexpansion'] [for pie plot]
        Labels for the categories defined by `groups`.
        
    ax : Optional[plt.Axes], optional
        Axes on which to plot. If None, a new figure will be created.
        
    figsize : Optional[Tuple[float, float]], optional
        Figure size, if `ax` is None.
        
    colors : Union[str, List[str], Dict[str, str]], optional, default='Set2'
        Color palette or list of colors.
        
    fontsize : int, optional, default=12
        Font size for labels and text.
        
    save : Optional[str], optional
        Path to save the figure. If None, the figure will not be saved.
        
    returnax : bool, optional, default=True
        Whether to return the Axes object.

    Returns
    -------
    Optional[plt.Axes]
        Returns the Axes object if `returnax` is True.
    """
#     adata.var['Chain'].unique()
    valid_list = ['TRA','TRB','IGH','IGK','IGL','IGHA','IGHG','IGHM','IGHD','IGHE']
    LiteralCheck(chain,valid_list)
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (3, 3)
        fig, ax = plt.subplots(figsize=figsize)


    if fontsize is None:
        fontsize = adjust_font_size(ax)

    if kind == 'bar':
        colors = get_colors(ntop, colors)
        expansion_df = expansion_tl(adata, chain=chain, focus_cdr3=focus_cdr3, groups=groups, labels=labels, show_match_dict=False)
        topdf = expansion_df.head(ntop)
        
        bars = ax.barh(topdf['ID'], topdf['n_cells'], edgecolor='black', color=colors,**kwarg)
        if label_loc != None:
            ax.bar_label(bars, labels=topdf['n_cells'], label_type=label_loc,padding = 0.2, fontsize=fontsize)
        ax.invert_yaxis()

    elif kind == 'pie':
        colors = get_colors(labels, colors)
        expansion_df = expansion_tl(adata, chain=chain, focus_cdr3=focus_cdr3, groups=groups, labels=labels, show_match_dict=False)
        statecount = expansion_df['group'].value_counts()
        statecount = statecount.loc[labels]
        total = expansion_df.shape[0]
        groups = statecount.index
        group_counts = statecount.values
        if 'edgecolor' not in kwarg.keys():
            edgecolor = 'black'
        else:
            edgecolor = kwarg['edgecolor']
            del kwarg['edgecolor']
        wedges = ax.pie(
            group_counts,
            startangle=90,
            radius=1,
            wedgeprops=dict(width=0.6, edgecolor=edgecolor),
            colors=colors,
            autopct=lambda p: f'{int(round(p * total / 100.0))}\n{p:.1f}%', 
            textprops={'fontsize': fontsize - 1, 'color': 'black'},
            pctdistance=0.7,
            **kwarg
        )
        ax.text(0, 0, f"{chain}\n{total}", ha='center', va='center', fontsize=fontsize)
        ax.legend(wedges[0], groups, bbox_to_anchor=(1, 1), frameon=False, fontsize=fontsize)

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')

    if returnax:
        return ax
    else:
        plt.show()


def hull_detect(
    adata: AnnData,
    inner: Union[List[str], str], 
    outter: Union[List[str], str],
    groupby: str,
    basis: str = 'spatial',
    cutoff: Optional[float] = None,
    ax: Optional[plt.Axes] = None, 
    figsize: Optional[Tuple[float, float]] = None,
    spotsize: float = 200,
    scatter_spotsize: Optional[float] = None,
    returnhull: bool = False,
    save: Optional[str] = None,
    base_alpha: float = 1,
    min_distance: Optional[float] = None,
    **kwarg,
) -> Optional[Union[plt.Axes, pd.DataFrame]]:
    """
    Detects the outer border cells based on proximity between inner and outer clusters in a spatial embedding.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    
    inner : Union[List[str], str]
        The inner clusters or groups.
    
    outter : Union[List[str], str]
        The outer clusters or groups.
    
    groupby : str
        The column in `adata.obs` to group the data by (e.g., clusters or cell types).
    
    basis : str, optional
        The embedding basis to use for spatial coordinates. Default is 'spatial'.
    
    cutoff : float, optional
        The percentile cutoff for the distance threshold. Default is 30.
    
    returnhull : bool, optional
        If True, returns the hull (border cells). If False, saves it to `adata.uns`. Default is False.

    spotsize : float, optional, default=200
        Spot size for base plot.

    scatter_spotsize : Optional[float], optional
        Spot size for scatter points.

        figsize : Optional[Tuple[float, float]], optional
        Figure size for the plot.

    save : Optional[str], optional
        File path to save the plot. If None, the plot is not saved.

    ax : Optional[plt.Axes], optional
        Axis to plot on. If None, a new figure is created.
    
    **kwarg
        kwarg for sns.scatterplot

    Returns
    -------
    pd.DataFrame or None
        Returns the hull (border cells) if `returnhull=True`, otherwise None (saves hull to `adata.uns`).
    """
    
    valid_obs = adata.obs.columns
    LiteralCheck(groupby, valid_obs)


    valid_cluster = adata.obs[groupby].unique().tolist()
    for ele in [inner, outter]:
        LiteralCheck(ele, valid_cluster)
    
    if isinstance(inner, str):
        inner = [inner]

    if isinstance(outter, str):
        outter = [outter]

    inner_site = inner
    set_inner = set(inner)
    set_outter = set(outter)
    set_intersection = set_inner.intersection(set_outter)
    if len(set_intersection) >= 1:
        print(f'Warning: Found intersection {list(set_intersection)}')
    
    adata.uns['immu_repertoire']['args']['hull'] = True
    hull_dict = {'inner': inner, 'outter': outter, 'groupby': groupby}

    # Extract spatial data for inner and outer cells
    inner_adata = adata[adata.obs[groupby].isin(inner)]
    outsite_adata = adata[adata.obs[groupby].isin(outter)] 
    inner = inner_adata.obsm[basis]
    outsite = outsite_adata.obsm[basis]
    
    # Build KDTree for inner group and query distances for outer group
    tree = KDTree(inner)
    result = tree.query(outsite, k=1)
    if cutoff == None and min_distance == None:
        min_distance = get_min_distance(adata,basis = basis)

    outsite_list = []
    for idx, dist in enumerate(result[0]):
        if cutoff != None:
            if dist[0] <= np.percentile(result[0], cutoff):
                outsite_list.append(idx)
        else:
            if dist[0] <= min_distance + 10 and dist[0] >= min_distance - 10:
                outsite_list.append(idx)
        
    border_cells = outsite_adata.obs.iloc[outsite_list].index
    hull = adata[border_cells].obsm[basis].copy()
    
    if 'immu_repertoire' not in adata.uns:
        adata.uns['immu_repertoire'] = {'args': {}}
    hull_dict['hull'] = hull
    adata.uns['immu_repertoire']['hull'] = hull_dict

    tree = KDTree(hull)
    obsdf = pd.DataFrame(adata.obsm[basis],index = adata.obs_names, columns = ['x','y'])
    obsdf[groupby] = adata.obs[groupby]
    distance,idxs = tree.query(obsdf[['x','y']],k = 1)
    obsdf['distance'] = distance
    obsdf['sign'] = obsdf[groupby].map(lambda x:  -1 if x in inner_site else 1)
    obsdf['abs_distance'] = obsdf['sign'].astype(int) * obsdf['distance']
    adata.obs['distance'] = obsdf['abs_distance']

    returnax = True
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (3, 3)
        fig, ax = plt.subplots(figsize=figsize)

    

    if scatter_spotsize == None:
        scatter_spotsize = spotsize/15
    sc.pl.spatial(adata,color = groupby,ax = ax,show =False,img_key = None,spot_size = spotsize, title = 'hull', basis = basis,alpha = base_alpha)
    sns.scatterplot(x = hull[:,0],y = hull[:,1],edgecolor = None,s = scatter_spotsize, **kwarg)
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if returnax:
        return ax
    
    fig.show()

def get_min_distance(adata,basis = 'spatial'):
    locations = adata.obsm[basis].copy()
    distances = np.linalg.norm(locations[:, np.newaxis] - locations, axis=2)
    np.fill_diagonal(distances, np.inf)
    min_distances = np.min(distances, axis=1)
    average_min_distance = np.mean(min_distances)
    return average_min_distance

def rank_distance(
    adata: AnnData,
    focus_cdr3: bool = False,
    clone_type: Optional[str] = None,
    min_spot: int = 10,
    ntop: int = 10,
    basis: str = 'spatial',
    colors: Dict[str, str] = {'top': '#fc5185', 'bottom': '#3fc1c9'},
    scatter_alpha: float = 0.2,
    scatter_color: str = 'grey',
    ax: Optional[plt.Axes] = None,
    figsize: Optional[Union[Tuple[float, float], None]] = None,
    save: Optional[str] = None,
    **kwarg,
) -> Optional[plt.Axes]:
    """
    Rank distances of cdr3s based on proximity to inner and outer cells and plot them.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    
    focus_cdr3 : bool, optional
        If True, focus on CDR3-based distances, otherwise on general clones.
    
    min_spot : int, optional
        Minimum number of cells (spots) required for distance calculation.
    
    ntop : int, optional
        Number of top and bottom elements to display.
    
    basis : str, optional
        The embedding basis to use for spatial coordinates. Default is 'spatial'.
    
    colors : dict, optional
        Dictionary specifying colors for 'top' and 'bottom' elements.
    
    scatter_alpha : float, optional
        Alpha transparency for the scatter plot points.
    
    scatter_color : str, optional
        Color for the scatter plot points.
    
    ax : plt.Axes, optional
        Matplotlib axis on which to plot. If None, a new figure and axis are created.
    
    figsize : Tuple[float, float], optional
        Size of the figure if `ax` is None.
    
    save : str, optional
        Path to save the figure. If None, the figure is not saved.
    
    **kwarg
        kwarg for sns.barplot
    Returns
    -------
    plt.Axes or None
        Returns the axis if `ax` is provided, otherwise returns None.
    """
    
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (3, 3)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        returnax = True

    if not adata.uns['immu_repertoire']['args'].get('hull', False):
        raise KeyError("Not found 'hull' in adata.uns, Please run pl.hull_detect first.")
    
    inner = adata.uns['immu_repertoire']['hull']['inner']
    outter = adata.uns['immu_repertoire']['hull']['outter']
    groupby = adata.uns['immu_repertoire']['hull']['groupby']

    inner_adata = adata[adata.obs[groupby].isin(inner)]

    tree = KDTree(adata.uns['immu_repertoire']['hull']['hull'])

    var_df = extract_var(adata, chain = clone_type, focus_cdr3 = focus_cdr3)

    if not focus_cdr3:
        clonedict = query(adata, var_names = var_df['ID'].tolist(), flatten=False)  
    else:
        clonedict = query(adata, cdr3_names = var_df['ID'].tolist(), flatten=False) 
        
    cdr3_total_distance = []
    for cdr3,cells in clonedict.items():
        xy = adata[cells].obsm[basis]
        if xy.shape[0] >= min_spot:
            distance_list = tree.query(xy)[0]
            sign_list = []
            for cell in cells:
                if cell in inner_adata.obs_names:
                    sign_list.append(-1)
                else:
                    sign_list.append(1)
            final_distance = []
            for distance,sign in zip(distance_list, sign_list):
                distance = distance[0] * sign
                final_distance.append(distance)
            for tmpd in final_distance:
                cdr3_total_distance.append([cdr3,tmpd])
    
    distancedf = pd.DataFrame(cdr3_total_distance, columns=['id', 'Distance'])

    adata.uns['immu_repertoire']['hull']['distance'] = distancedf
    print("Save distance in .uns['immu_repertoire']['hull']['distance']")

    median_distance = distancedf.groupby('id')['Distance'].mean()
    sorted_ids = median_distance.sort_values().index
    top_ids = sorted_ids[:ntop]
    bottom_ids = sorted_ids[-ntop:]
    
    selected_ids = pd.Index(list(top_ids) + list(bottom_ids))
    filtered_df = distancedf[distancedf['id'].isin(selected_ids)].copy()
    
    palette = {id_: 'top' for id_ in top_ids}
    palette.update({id_: 'bottom' for id_ in bottom_ids})
    filtered_df['color'] = filtered_df['id'].map(palette)
    if focus_cdr3:
        rename = 'CDR3'
    else:
        rename = 'Clone'
    filtered_df = filtered_df.rename({'id':rename},axis = 1)
    sns.barplot(data=filtered_df, x='Distance', y=rename, order=selected_ids, palette={'top': colors['top'], 'bottom': colors['bottom']}, hue='color', ax=ax, **kwarg)
    sns.stripplot(data=filtered_df, x='Distance', y=rename, color=scatter_color, alpha=scatter_alpha, ax=ax)
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if returnax:
        return ax
    else:
        plt.show()