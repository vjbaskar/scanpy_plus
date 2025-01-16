from ..globimport import *

def histplot(data, keys, groupby=None, xaxis_all_visibility = True, **kwargs):
    import scanpy as sc
    """
    Similar to scanpy violin plot but here histograms are plotted
    
    Params:
    data: anndata obj
    keys: str. obs column that should be plotted
    groupby: facet grid column to be used
    xaxis_all_visibility: By default X axis tick labels are visible in bottom-most plots. Make it visible in all plots instead
    kwargs: Args to be passed to sc.pl.umap
    
    Returns: None
    
    Note: You have to have colours for the obs col defined in anndata.uns
    
    Example usage:
    cf.pl.histplot(adata, keys='log1p_total_counts', groupby='samplename', col_wrap = 4)
    
    
    """
    #data = rna
    df = data.obs
    #keys = 'log1p_total_counts'
    #groupby = 'condition'
    if groupby == None:
        g = sns.displot(df, x = keys)
    else:
        g = sns.displot(df, x = keys, col = groupby, facet_kws=dict(sharey=False), hue=groupby, **kwargs)
    for ax in g.axes:
        ax.tick_params(labelbottom=xaxis_all_visibility)
    return g