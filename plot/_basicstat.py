import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData
from typing import Union, List, Optional, Dict, Literal
import numpy as np
import pandas as pd 

from venn import venn
from upsetplot import UpSet,from_contents


from .._utils import valid_chain, LiteralCheck, extract_var, GroupbyCheck
from ._utils import get_colors,adjust_font_size
from ..tools import query
from ..process import filter_clone
##############################
## CDR3 length
##############################

def cdr3_length_histogram(
    adata: AnnData,
    chains: Union[str, List[str]] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Optional[tuple] = None,
    colors: Optional[Union[str, List[str], Dict[str, str]]] = None,
    alpha: float = 1.0,
    legend_loc: Optional[tuple] = (1, 1),
    save: Optional[str] = None,
    **kwag
):
    """
    Plot histogram of CDR3 lengths for specified chains in the given AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing clone information.

    chains : Union[str, List[str], None], optional, default=None
        Chains to plot. Can be a string or a list of strings. If None, all valid chains are plotted.

    ax : matplotlib.axes.Axes, optional, default=None
        The axis to plot on. Only used if `show` is False.

    figsize : tuple, optional, default=None
        Figure size for the plot, used only if `show` is True and `ax` is None.
        If None, defaults to (5, 5).
    
    colors : Union[str, List[str], Dict, None], optional, default=None
        If colors is None, the default is tab10. If colors is str, the color is obtained through plt.get_cmap().
        If colors is a list, it needs to match the number of elements to be drawn.

    alpha : float, optional, default=1
        Alpha transparency level for the histogram bars.

    legend_loc : tuple, optional, default=(1, 1)
        The location of the legend on the plot. If None, no legend is displayed.

    save : str, optional, default=None
        File path to save the plot. If None, the plot is not saved.

    **kwag : dict, optional
        Additional keyword arguments for `plt.hist`

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axis containing the plot, returned if `show=False` and `save=None`.
    """

    # Chain selection
    returnax = True
    if chains is None:
        chains = valid_chain(adata, isotype=False)
    else:
        total_chains = valid_chain(adata, isotype=True)
        LiteralCheck(chains, total_chains)
    
    colors = get_colors(chains,colors)

    if 'edgecolor' not in kwag:
        kwag['edgecolor'] = 'black'
    
    adata = filter_clone(adata,min_cells = 1,min_cdr3_length = 0 ,max_cdr3_length = 100000)
    clone_var_df = adata.var.copy()

    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (4, 4)
        fig, ax = plt.subplots(figsize = figsize)
    

    for chain,color in zip(chains,colors):
        if len(chain) == 3:
            subclone_var = clone_var_df[clone_var_df['Chain'] == chain]
        elif len(chain) == 4:
            subclone_var = clone_var_df[clone_var_df['IGHisotype'] == chain]
        
        cdr3_lengths = subclone_var['CDR3_lenth']
        ax.hist(cdr3_lengths, bins=range(min(cdr3_lengths), max(cdr3_lengths) + 1), alpha=alpha, label=chain, color = color, **kwag)

    if legend_loc is not None:
        ax.legend(bbox_to_anchor=legend_loc)

    if save != None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    if returnax:
        return ax
    fig.show()

##############################
## Clone stat or Clone doninate gene
##############################

def gene_distribution(
    adata: AnnData,
    region: Literal['V', 'D', 'J', 'C'],
    ntop: int = 5,
    order: Optional[Union[str, List[str]]] = None,
    colors: Optional[Union[str, List[str], Dict[str, str]]] = None,
    othercolor: str = 'grey',
    ax: Optional[plt.Axes] = None,
    figsize: Optional[tuple] = None,
    save: Optional[str] = None,
    returndf: bool = False,
    **kwarg
):
    """
    Plot a stacked bar plot of region hits (e.g., V, D, J, C) for clone chains.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing clone information.

    region : str, optional, default='V'
        The region to analyze and plot. Can be one of ['V', 'D', 'J', 'C'].

    ntop : int, optional, default=5
        The number of top genes to highlight in the plot.

    order : Optional[Union[str, List[str]]], optional, default=None
        The order of the abscissa of the histogram, used to sort or filter different chains. If None, no order is applied.

    colors : Optional[Union[str, List[str], Dict]], optional, default=None
        Colors for the top ntop genes. Can be a single color, a list of colors, or a dictionary mapping genes to colors.

    othercolor : str, optional, default='grey'
        Color for genes not in the top ntop list.

    ax : Optional[plt.Axes], optional, default=None
        The matplotlib axis to plot on. If None, a new figure and axis are created.

    figsize : Optional[tuple], optional, default=None
        The size of the figure. Used if ax is None and show is True.

    save : Optional[str], optional, default=None
        The file path to save the plot. If None, the plot is not saved.

    **kwarg : dict
        Additional keyword arguments for the bar plot (e.g., edgecolor, linewidth).

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axis containing the plot, returned if `save=None`.

    """
    clone_df = adata.uns['immu_repertoire']['clone_info']
    clone_df = clone_df.replace('None',np.nan)
    LiteralCheck(region, ['V', 'D', 'J', 'C'])
    if f'{region}_align' not in clone_df.columns:
        raise KeyError(f'Clone info not contain {region}_align')

    result_list = []
    for chain, tmpdf in clone_df.groupby('topChains'):
        vc = tmpdf[f'{region}_align'].value_counts()
        result_list.extend(list(zip(vc.index, vc.values, [chain] * vc.shape[0])))

    df = pd.DataFrame(result_list, columns=[f'{region}Hits', 'count', 'chain'])
    
    df = df[df[f'{region}Hits'] != 'NA']
    df['rank'] = df.groupby('chain')['count'].rank(method='first', ascending=False)

    if not isinstance(colors,dict):
        colors = get_colors(['1'] * ntop, colors=colors)
        df['color'] = df.apply(lambda row: colors[int(row['rank']) - 1] if row['rank'] <= ntop else othercolor, axis=1)
    else:
        df['color'] = df.apply(lambda row: colors[row[f'{region}Hits']] if row['rank'] <= ntop else othercolor, axis=1)
    if returndf:
            return df
    returnax = True
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (4, 4)
        fig, ax = plt.subplots(figsize = figsize)

    chains = df['chain'].unique()
    bardict = dict(df.groupby('chain')[f'{region}Hits'].agg(lambda x: len(list(x))))

    if order is not None:
        LiteralCheck(order, chains)
        chains = [i for i in order if i in chains]

    for chain in chains:
        chain_data = df[df['chain'] == chain].sort_values(by='count', ascending=False)
        bottom = None
        for _, row in chain_data.iterrows():
            ax.bar(chain, row['count'], bottom=bottom, color=row['color'], 
                label=f'Rank {int(row["rank"])}' if row['rank'] <= ntop else "", **kwarg)
            bottom = row['count'] if bottom is None else bottom + row['count']

    new_labels = [f'{chain}:{bardict[chain]}' if chain in bardict else chain for chain in chains]
    ax.set_xticks(range(len(chains)))
    ax.set_xticklabels(new_labels)
    ax.set_ylabel('Clone Nums')

    handles, labels = [], []
    for chain in chains:
        chain_data = df[df['chain'] == chain]
        top_5_chain_data = chain_data[chain_data['rank'] <= ntop]  # Only top 5
        for _, row in top_5_chain_data.iterrows():
            handles.append(plt.Rectangle((0, 0), 1, 1, color=row['color']))
            labels.append(f"{row[f'{region}Hits']} ({chain})")

    ncol = 1 if min(bardict.values()) < ntop else len(bardict.keys())
    ax.legend(handles, labels, title=f'{region} Gene by Chain', bbox_to_anchor=(1.05, 1), loc='upper left', ncol=ncol)


    ax.grid(False)
    if save != None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    if returnax:
        return ax
    fig.show()

##############################
## Clone stat or Clone doninate gene
##############################

def _absolute_value(pct,total):
    return int(round(pct*total/100.0))

def clone_distribution(
    adata: AnnData,
    isotype: bool = False,
    focus_cdr3: bool = False,
    order: Optional[Union[List[str], None]] = None,
    kind: Literal['pie', 'bar'] = 'bar',
    colors: Optional[Union[str, List[str], Dict, None]] = None,
    ax: Optional[plt.Axes] = None,
    figsize: Optional[tuple] = None,
    save: Optional[str] = None,
    fontsize: Optional[int] = None,
    pie_width: float = 0.5,
    **kwarg,
):
    """
    Plot the distribution of BCR or TCR chains using either a bar plot or pie chart.

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object containing clone information.
    isotype : bool, optional
        Whether to consider BCR isotypes (default: False).
    focus_cdr3 : bool, optional
        Whether to focus on CDR3 sequence, not clone (default: False).
    order : Optional[List[str], None], optional
        The order of the abscissa of the histogram, used to sort or filter different chains (default: None).
    kind : {'pie', 'bar'}, optional
        Type of plot to display (default: 'bar').
    colors : Optional[Union[str, List[str], Dict, None]], optional
        Colors for the chains (default: None).
    ax : Optional[plt.Axes], optional
        Existing axes to plot on (default: None).
    figsize : Optional[tuple], optional
        Size of the figure (default: None).
    save : Optional[str], optional
        Path to save the figure (default: None).
    fontsize : int, optional
        Font size for labels (default: None).
    pie_width : float, optional
        Width of pie chart wedges (default: 0.5).
    **kwarg : dict
        Additional keyword arguments for `sns.barplot`.

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object with the plot if `ax` is None, otherwise returns `None`.
    """
    
    returnax = True
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (5, 5)
        fig, ax = plt.subplots(figsize=figsize)

    var_df = extract_var(adata, focus_cdr3 = focus_cdr3)
    
    if not isinstance(isotype,bool):
        raise KeyError(f'isotype must be bool')
    if not isotype:
        chain_counts = var_df['Chain'].value_counts()
    else:
        if 'IGHisotype' not in var_df.columns:
            raise KeyError('Not Found IGHisotype in adata.var. Please enrich your clone information')
        var_df['IGHisotype'] = var_df['IGHisotype'].replace('NoBCR_C',np.nan)
        chain_counts = var_df['IGHisotype'].value_counts()

    chains_for_plot = chain_counts.index.tolist()
    colors = get_colors(chains_for_plot, colors=colors)

    if fontsize == None:
        fontsize = adjust_font_size(ax)
    
    if kind == 'bar':
        barplot = sns.barplot(x = chain_counts.index, y = chain_counts.values, hue = chain_counts.index, palette=colors, edgecolor='black',order = order, ax=ax, **kwarg)
        for p in barplot.patches:
            barplot.annotate(format(p.get_height(), '.0f'), 
                             (p.get_x() + p.get_width() / 2., p.get_height()), 
                             ha='center', va='center', xytext=(0, 5), fontsize = fontsize, 
                             textcoords='offset points')

    elif kind == 'pie':
        total = chain_counts.sum()
        plot_value = chain_counts.values.tolist()
        wedges = ax.pie(
            plot_value,
            startangle=90,
            radius=1,
            wedgeprops=dict(width=pie_width, edgecolor='black'),
            colors=colors,
            autopct=lambda p: f'{_absolute_value(p, total)}',
            textprops={'fontsize': fontsize, 'color': 'black'},
            pctdistance=(1 - pie_width) / 2 + pie_width,
        )
        tmpname = 'Clone'
        if focus_cdr3:
            tmpname = 'CDR3'
        ax.text(0, 0, f"{tmpname}\n{total}", ha='center', va='center', fontsize=fontsize)
        ax.legend(wedges[0], chain_counts.index, loc='upper left', bbox_to_anchor=(0.95, 0.95), frameon=False)
    else:
        raise KeyError("Kind must be in ['pie','bar']")
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')

    if returnax:
        return ax
    else:
        plt.show()

########################################
#####  intersection
########################################

def intersection(
    adata: AnnData, 
    groupby: str, 
    chain: Optional[str] = None, 
    focus_cdr3: bool = False, 
    groups: Optional[Union[str, List[str]]] = None, 
    kind: Literal['venn', 'upset'] = 'venn',
    ax: Optional[plt.Axes] = None,
    figsize: Optional[tuple] = None,
    save: Optional[str] = None,
    **kwarg
):

    returnax = True
    if kind == 'venn':
        if ax is None:
            returnax = False
            if figsize is None:
                figsize = (5, 5)
            fig, ax = plt.subplots(figsize=figsize)

    candidate = GroupbyCheck(adata, groupby = groupby, groups = groups)
    
    cdr3_dict = {}
    
    for i in candidate:
        cells = adata[adata.obs[groupby] == i].obs_names
        region_cdr3s = query(adata, obs_names=cells, obs_chain=chain, focus_cdr3=focus_cdr3, flatten=True)
        cdr3_dict[i] = set(region_cdr3s)
    
    if kind == 'venn':
        g = venn(cdr3_dict,ax = ax,**kwarg)
    elif kind == 'upset':
        upset_nowarming(cdr3_dict,**kwarg)
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')

    if kind == 'venn' and returnax:
        return ax
    else:
        plt.show()
    
def upset_nowarming(cdr3_dict, **kwarg):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        inputupset = from_contents(cdr3_dict)
        UpSet(inputupset, subset_size="count",**kwarg).plot()

