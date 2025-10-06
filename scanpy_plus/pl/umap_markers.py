from ..globimport import *

def umap_markers(adata, markers, **kwargs):
    """
    UMAP of a set of markers
    adata: anndata
    markers: a dict of markers. pd.DataFrame({"celltype1":[Gene1,Gene2...]})
    kwargs: passed to sc.pl.umap
    """
    for k,v in markers.items():
        olap = list(set(adata.var_names) & set(v))
        if len(olap) > 0:
            sc.pl.umap(adata, color=olap, ncols=len(olap), title = list(map(lambda x: "/".join([k,x]), olap)), **kwargs)






