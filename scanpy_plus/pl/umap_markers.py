from ..globimport import *

def haem_markers(adata, markers, **kwargs):
    """
    """
    for k,v in markers.items():
        olap = list(set(adata.var_names) & set(v))
        if len(olap) > 0:
            sc.pl.umap(adata, color=olap, ncols=len(olap), title = list(map(lambda x: "/".join([k,x]), olap)), **kwargs)






