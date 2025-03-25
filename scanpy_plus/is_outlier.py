from .globimport import *
import numpy as np

def is_outlier(adata, 
               metric: str, 
               nmads: int) -> list:
    """
    Function from best practices for scrnaseq handbook
    adata: anndata obj
    metric: obs column (string or list)
    nmads: mads (int)
    
    return: list of cell (obs) indices
    """
    from scipy.stats import median_abs_deviation
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier