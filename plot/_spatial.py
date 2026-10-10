from typing import Union, List, Optional, Dict, Literal, Tuple
import matplotlib.pyplot as plt
from anndata import AnnData
import seaborn as sns
import numpy as np
import scanpy as sc 
import pandas as pd



from ..tools import query
from .._utils import LiteralCheck

from ._utils import get_colors, generate_gradient_colors


def spatial(
    adata: AnnData,
    clones: Optional[Union[str, List[str], Dict[str, List[str]]]] = None,
    cdr3s: Optional[Union[str, List[str], Dict[str, List[str]]]] = None,
    obs: Optional[Union[str, List[str]]] = None,
    base: Optional[Literal['image', str]] = None,
    spotsize: float = 200,
    scatter_spotsize: Optional[float] = None,
    no_scatter: bool =  False,
    colors: Optional[Union[str, List[str], Dict[str, str]]] = None,
    basis: str = 'spatial',
    base_alpha: float = 0.5,
    scatter_alpha: float = 1.0,
    kde: bool = False,
    fill: bool = True,
    kde_bw_adjust: Union[List[float], float] = 0.5,
    kde_thresh: Union[List[float], float] = 0.5,
    kde_alpha: Union[List[float], float] = 0.5,
    kde_fill: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
    save: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    legend_loc: Optional[Tuple[float, float]] = (1, 1),
    light_factor: Optional[float] = None,
    na_color: Optional[str] = None,
) -> Optional[plt.Axes]:
    """
    Plot clones or CDR3s spatial distribution on a scatter plot, with optional KDE.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing clone or CDR3 information.
    
    clones : Optional[Union[str, List[str], Dict[str, List[str]]]], optional
        Clone names to plot. Can be a string, a list of strings, or a dictionary where keys are clone names
        and values are lists of sequences.

    cdr3s : Optional[Union[str, List[str], Dict[str, List[str]]]], optional
        CDR3 sequences to plot. Can be a string, a list of strings, or a dictionary where keys are CDR3 names
        and values are lists of sequences.

    base : Literal['image', str], optional, default='image'
        A column in `adata.obs` or 'image' for a background.

    spotsize : float, optional, default=200
        Spot size for base plot.

    scatter_spotsize : Optional[float], optional
        Spot size for scatter points.

    colors : Optional[Union[str, List[str], Dict[str, str]]], optional
        Colors for the different clones or CDR3s.

    basis : str, optional, default='spatial'
        Basis for spatial plotting, e.g., 'spatial' coordinates.

    base_alpha : float, optional, default=0.5
        Alpha transparency for the base image or plot.

    scatter_alpha : float, optional, default=1.0
        Alpha transparency for scatter points.

    kde : bool, optional, default=False
        Whether to plot a KDE over the scatter points.

    kde_bw_adjust : float, optional, default=0.5
        KDE bandwidth adjustment.

    fill : bool, optional, default=True
        Whether to fill the KDE plot.

    kde_thresh : float, optional, default=0.5
        KDE threshold for plotting.

    kde_alpha : float, optional, default=0.5
        Alpha transparency for the KDE plot.

    kde_fill : bool, optional, default=True
        Whether to fill the KDE area.

    figsize : Optional[Tuple[float, float]], optional
        Figure size for the plot.

    save : Optional[str], optional
        File path to save the plot. If None, the plot is not saved.

    ax : Optional[plt.Axes], optional
        Axis to plot on. If None, a new figure is created.

    legend_loc : Optional[Tuple[float, float]], optional, default=(1, 1)
        Location of the legend.

    Returns
    -------
    Optional[plt.Axes]
        The axis containing the plot, returned if `ax` is not None.
    """
    
    # Ensure either clones or CDR3s is given, not both
    if clones is not None and cdr3s is not None:
        raise KeyError('Both clones and cdr3s were given.')

    returnax = True
    if ax is None:
        returnax = False
        if figsize is None:
            figsize = (4, 4)
        fig, ax = plt.subplots(figsize=figsize)

    scalefactor = baseplot(adata, ax=ax, base=base, base_alpha=base_alpha, spotsize=spotsize, basis=basis,na_color = na_color, )

    if scatter_spotsize == None:
        scatter_spotsize = spotsize/15

    if clones is not None:
        if not isinstance(clones, dict):
            clonesdict = {'Clone(s)': clones}
        else:
            clonesdict = clones
        

        colors = get_colors(clonesdict.keys(), colors)

        cells_list = []
        for (clonename, clones), color in zip(clonesdict.items(), colors):
            cells = query(adata, var_names=clones, flatten=True)
            cells_list.append(cells)
            if not no_scatter:
                scatterplot(adata, cells=cells, scalefactor=scalefactor, ax=ax, scatter_spotsize=scatter_spotsize, scatter_alpha=scatter_alpha, label=clonename, color=color, basis = basis)

        if kde:
            n = len(clonesdict)
            kde_bw_adjusts,kde_threshs,kde_alphas = adjust_kde_argv(kde_bw_adjust,kde_thresh,kde_alpha,n)
            for (clonename, clones), color, cells, kde_bw_adjust, kde_thresh, kde_alpha in zip(clonesdict.items(), colors, cells_list,kde_bw_adjusts,kde_threshs,kde_alphas):
                kdeplot(adata, cells=cells, scalefactor=scalefactor, ax=ax, bw_adjust=kde_bw_adjust, fill=kde_fill, fill_color=color, kde_thresh=kde_thresh, kde_alpha=kde_alpha, basis = basis, light_factor = light_factor)

    elif cdr3s is not None:
        if not isinstance(cdr3s, dict):
            cdr3sdict = {'CDR3(s)': cdr3s}
        else:
            cdr3sdict = cdr3s

        colors = get_colors(cdr3sdict.keys(), colors)

        cells_list = []
        for (cdr3name, cdr3s), color in zip(cdr3sdict.items(), colors):
            cells = query(adata, cdr3_names=cdr3s, flatten=True)
            cells_list.append(cells)
            if not no_scatter:
                scatterplot(adata, cells=cells, scalefactor=scalefactor, ax=ax, scatter_spotsize=scatter_spotsize, scatter_alpha=scatter_alpha, label=cdr3name, color=color, basis = basis)

        if kde:
            n = len(cdr3sdict)
            kde_bw_adjusts,kde_threshs,kde_alphas = adjust_kde_argv(kde_bw_adjust,kde_thresh,kde_alpha,n)
            for (cdr3name, clones), color, cells, kde_bw_adjust, kde_thresh, kde_alpha in zip(cdr3sdict.items(), colors, cells_list,kde_bw_adjusts,kde_threshs,kde_alphas):
                kdeplot(adata, cells=cells, scalefactor=scalefactor, ax=ax, bw_adjust=kde_bw_adjust, fill=kde_fill, fill_color=color, kde_thresh=kde_thresh, kde_alpha=kde_alpha, basis = basis, light_factor = light_factor)
    elif obs is not None:
        # when obs  only plot kde
        if not isinstance(obs, list):
            obs = [obs]
        
        colors = get_colors(obs, colors)
        print("When specific obs, only plot KDE, if need spatial plot, please use 'sc.pl.spatial' ")
        n = len(obs)
        kde_bw_adjusts,kde_threshs,kde_alphas = adjust_kde_argv(kde_bw_adjust,kde_thresh,kde_alpha,n)
        for obs, color, kde_bw_adjust, kde_thresh, kde_alpha in zip(obs, colors, kde_bw_adjusts,kde_threshs,kde_alphas):
            kdeplot(adata, value_obs = obs, scalefactor=scalefactor, ax=ax, bw_adjust=kde_bw_adjust, fill=kde_fill, fill_color=color, kde_thresh=kde_thresh, kde_alpha=kde_alpha, basis = basis, light_factor = light_factor)

    
    if legend_loc is not None:
        ax.legend(bbox_to_anchor=legend_loc, loc='best', frameon=False)

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    
    if returnax:
        return ax
    
    fig.show()

def baseplot(
    adata: AnnData,
    ax: plt.Axes,
    base: Optional[Literal['image',str]] = None, 
    base_alpha: float = 1,
    spotsize: float = 200,
    basis: str = 'spatial',    
    na_color: str = "#D3D3D3",
):
    data = adata.copy()
    scalefactor = 1
    if base == 'image':
        if 'spatial' not in data.uns.keys():
            raise KeyError("Not found image infomation in adata.uns['spatial']")
        libname = list(data.uns['spatial'].keys())[0]
        image = data.uns['spatial'][libname]['images']['hires']
        scalefactor = data.uns['spatial'][libname]['scalefactors']['tissue_hires_scalef']
        data.obs['_tmp_plot_always_not_found_'] = np.nan
        tmpfig, tmpax = plt.subplots(figsize=(1, 1))
        sc.pl.spatial(data, color = '_tmp_plot_always_not_found_', spot_size = 1, show=False, ax = tmpax)
        spatial_xlim = tmpax.get_xlim()
        spatial_ylim = tmpax.get_ylim()
        plt.close()
        ax.imshow(image, alpha = base_alpha)
        ax.set_xlim(spatial_xlim)
        ax.set_ylim(spatial_ylim)
        ax.grid(False)
    elif base == None:
        data.obs[f'loci'] = 'NA'
        if na_color is None:
            na_color = '#D3D3D3'
        data.uns[f'loci_colors'] = [na_color]
        sc.pl.spatial(data, color = 'loci', spot_size = spotsize, basis = basis, show = False, ax = ax, legend_loc = None,img_key = None, alpha = base_alpha, title=None)
    else:
        LiteralCheck(base,list(data.obs.columns))
        sc.pl.spatial(data, color = base, spot_size = spotsize, basis = basis, legend_loc = None, show = False,img_key = None, ax = ax, alpha = base_alpha, title="")
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
    return scalefactor

def adjust_kde_argv(kde_bw_adjust,kde_thresh,kde_alpha,n):
    if not isinstance(kde_bw_adjust, list):
        kde_bw_adjust = [kde_bw_adjust]*n
    if not isinstance(kde_thresh, list):
        kde_thresh = [kde_thresh]*n
    if not isinstance(kde_alpha, list):
        kde_alpha = [kde_alpha]*n
    return kde_bw_adjust,kde_thresh,kde_alpha


def scatterplot(
    adata: AnnData,
    cells: List[str], 
    scalefactor: float, 
    ax: plt.Axes,
    scatter_spotsize: float = 200,
    label: Optional[str] = None,
    color: Optional[str] = None,
    scatter_alpha: float = 0.5,
    basis: str = 'spatial',
    **kwarg
):
    basis_df = pd.DataFrame(adata[cells].obsm[basis],columns = ['x', 'y'])
    basis_df['x'] *= scalefactor
    basis_df['y'] *= scalefactor
    if 'edgecolor' not in kwarg.keys():
        kwarg['edgecolor'] = None
    ax.scatter(data = basis_df, x = 'x', y = 'y', s = scatter_spotsize, alpha = scatter_alpha, label = label, color = color, **kwarg)
    

def kdeplot(
    adata: AnnData,
    scalefactor: float, 
    ax: plt.Axes,
    value_obs: Optional[str] = None,
    cells: Optional[List[str]] = None, 
    bw_adjust: float = 0.3,
    fill: bool = True,
    fill_color: Optional[str] = None,
    kde_thresh: float = 0.5,
    kde_alpha: float = 1,
    label: Optional[str] = None,
    basis: str = 'spatial',
    light_factor: Optional[float] = None,
):
    if cells != None:
        basis_df = pd.DataFrame(adata[cells].obsm[basis],columns = ['x', 'y'])
    else:
        basis_df = pd.DataFrame(adata.obsm[basis],columns = ['x', 'y'])

    
    basis_df['x'] *= scalefactor
    basis_df['y'] *= scalefactor
    

    # light color :
    kde_argv = {}
    if light_factor == None:
        kde_argv['color'] = fill_color
    else:
        fill_color = generate_gradient_colors(basecolor = fill_color, ncolor = None, light_factor = light_factor, reverse = True)
        kde_argv['cmap'] = fill_color

    if value_obs != None:
        basis_df['value'] = adata.obs[value_obs].tolist()
        if basis_df['value'].max() == basis_df['value'].min():
            pass
        else:
            basis_df['value'] = (basis_df['value'] - basis_df['value'].min()) / (basis_df['value'].max() - basis_df['value'].min())
        contour = sns.kdeplot(data = basis_df, x = 'x', y = 'y',weights = 'value', bw_adjust = bw_adjust, fill = fill, thresh = kde_thresh, alpha = kde_alpha, **kde_argv)

    else:
        contour = sns.kdeplot(data = basis_df, x = 'x', y = 'y', bw_adjust = bw_adjust, fill = fill, thresh = kde_thresh, alpha = kde_alpha, **kde_argv)

    fig = ax.get_figure()
    cbar = fig.colorbar(contour.collections[-1], ax=ax, label='Density', fraction=0.05, pad=0.1,aspect = 10)
    cbar.set_ticks([cbar.vmin, cbar.vmax])
    if value_obs ==None:
        max_label = "more\nclusterd"
        min_label = "less"
        cbar.set_ticklabels([min_label, max_label]) 
    else:
        pass
    cbar.set_label(label)