from .globimport import *


def rankobs(adata, obs = ['total_counts', 'n_genes_by_counts', 'doublet_score', 'pct_counts_mt'], ascending = [True, True, False, False], percentile = True):
    """
    Rank observations:
    
    adata: anndata
    obs: list of string. obs cols to be ranked
    ascending: bool. Ascending or descending order
    percentile: bool. compute in percs
    
    return: dataframe
    
    """
    obsdf = adata.obs
    df = pd.DataFrame()
                          
    for obsvar, high in zip(obs, ascending):
        df[obsvar] = obsdf[obsvar].rank(ascending = high, pct = percentile)
    return df
