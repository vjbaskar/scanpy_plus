from ..globimport import *

def splitplot(adata, clust_key, colorby, size=60, basis='X_umap', frameon=False, legend_loc=None, **kwargs):
    import scanpy as sc
    """
    Plot each factor of an obs column as separate plots.
    
    Params:
    adata: anndata obj
    clust_key: str. obs column that should be used for splitting
    size: int. Point size
    frameon: bool.
    legend_loc: str. Same as in sc.pl.umap
    kwargs: Args to be passed to sc.pl.umap
    
    Returns: None
    
    Note: You have to have colours for the obs col defined in anndata.uns
    
    
    
    """
    
    tmp = adata.copy()
    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]
    #sc.pl.embedding(tmp, basis = basis, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)
    sc.pl.embedding(tmp, basis = basis, groups=tmp.obs[clust].cat.categories[1:].values, color=colorby, size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

