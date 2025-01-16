from ..globimport import *

def splitscatter(data, x, y, groupby=None, colorby=None, **kwargs):
    import scanpy as sc
    """
    Similar to scanpy violin plot but here histograms are plotted
    
    Params:
    data: anndata obj
    keys: str. obs column that should be plotted
    groupby: facet grid column to be used
    kwargs: Args to be passed to sc.pl.umap
    
    Returns: None
    
    Note: You have to have colours for the obs col defined in anndata.uns
    
    
    
    """
    #data = rna
    df = data.obs
    #keys = 'log1p_total_counts'
    #groupby = 'condition'
    if groupby == None:
        g = sns.scatterplot(df, x = x, y = y, **kwargs)
    else:
        if colorby == None: colorby == groupby
        g = sns.displot(df, x = x, y = y, col = colorby, facet_kws=dict(sharey=False), hue=groupby, **kwargs)
    return g