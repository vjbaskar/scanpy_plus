from .globimport import *

def get_from_raw(data):
    """
    Does anndata.raw.to_adata() in a more streamlined way.
    The to_data() method does not carry all the metadata such as uns, varm etc
    This one faithfully keeps the genes the same and also adds in uns
    
    Requires: anndata.raw to be set
    
    Returns: anndata
    
    Example: get_from_raw(adata_scaled)
    
    """
    temp = data.raw.to_adata()
    x = data.copy()
    x.X = temp[:, data.var_names].X
    x.uns['log1p']['base'] = None
    return x