from ..globimport import *


def rankobs(adata, obs = ['total_counts', 'n_genes_by_counts', 'doublet_score', 'pct_counts_mt'], ascending = [True, True, False, False], percentile = True):
    """
    Rank observations.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.
    obs : list of str
        List of observation column names to be ranked.
    ascending : bool, optional
        Sort in ascending order (default is True).
    percentile : bool, optional
        If True, compute rankings in percentiles (default is False).

    Returns
    -------
    pd.DataFrame
        A DataFrame with ranked observation values.
    """
    
    
    obsdf = adata.obs
    df = pd.DataFrame()
                          
    for obsvar, high in zip(obs, ascending):
        df[obsvar] = obsdf[obsvar].rank(ascending = high, pct = percentile)
    return df
