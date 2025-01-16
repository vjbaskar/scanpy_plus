from ..globimport import *




def splitplot(adata, split_by, colour_by, size=60, frameon=False, legend_loc=None, **kwargs):
    import scanpy as sc
    """
    Plot each factor of an obs column as separate plots.
    
    Params:
    adata: anndata obj
    split_by: str. obs column that should be used for splitting
    colour_by: str. obs column that should be used for colouring within each split. Can be same as split_by
    size: int. Point size
    frameon: bool.
    legend_loc: str. Same as in sc.pl.umap
    kwargs: Args to be passed to sc.pl.umap
    
    Returns: None
    
    Note: You have to have colours for the obs col defined in anndata.uns
    
    
    
    """
    
    tmp = adata.copy()
   # for i,clust in enumerate(adata.obs[split_by].cat.categories):
   #     #print(clust)
   #     tmp.obs[clust] = adata.obs[split_by].isin([clust]).astype('category')
   #     tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[split_by+'_colors'][i]]
   # sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=tmp.obs[split_by].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)
  
    for i,clust in enumerate(adata.obs[split_by].cat.categories):
        tmp.obs[clust] = adata.obs[colour_by].astype(str)
        tmp.obs.loc[ ~ tmp.obs[split_by].isin([clust]), clust] = None
        tmp.obs[clust] = tmp.obs[clust].astype("category")
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[split_by+'_colors'][i]]
    print(tmp.obs.head())
    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=tmp.obs[split_by].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)
        
    return tmp

def splitplot(adata, split_by = 'condition', colour_by = 'celltype', **kwargs):
    
    """
    Plot each factor of an obs column as separate plots.
    
    Params:
    adata: anndata obj
    split_by: str. obs column that should be used for splitting
    colour_by: str. obs column that should be used for colouring within each split. Can be same as split_by
    size: int. Point size
    frameon: bool.
    legend_loc: str. Same as in sc.pl.umap
    kwargs: Args to be passed to sc.pl.umap
    
    Returns: None
    
    Note: You have to have colours for the obs col defined in anndata.uns
    
    
    
    """
    data = adata.copy()
    #split_by = 'condition'
    #colour_by = 'celltype'
    split_vals = data.obs[split_by].unique()
    #print(split_vals)
    for split in split_vals:
        #split = 'LK'
        #print(f"***{split}", end = "")
        data.obs[split] = data.obs[colour_by]
        data.obs.loc[~data.obs[split_by].isin([split]), split] = None

        data.obs[split] = data.obs[split].astype('category')
    
        data.uns[split+'_colors'] = ['#d3d3d3']  + data.uns[colour_by +'_colors']
        data.uns[split+'_colors'] = data.uns[colour_by +'_colors']
        
    #print(data.obs.head())
    sc.pl.umap(data,color=split_vals, **kwargs)
