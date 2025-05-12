from ..globimport import *
import scipy
import seaborn as sns
def assign_sex(data, ygenes, use_raw = True):
    """
    Assigns sex based on Xist and Y-chromosome gene expression.

    This function expects log-normalised input. It adds the following columns to `adata.obs`:
    
    - `xist_logn`: log-normalized expression of Xist
    - `Ygene_logn`: mean log-normalized expression of the specified Y-chromosome genes
    - `xist_bin`: binary indicator for Xist expression
    - `Ygene_bin`: binary indicator for Y gene expression
    - `sex`: inferred sex ("male", "female", or "ambiguous")

    Parameters
    ----------
    data : anndata.AnnData
        Annotated data matrix.
    ygenes : list of str
        List of Y-chromosome gene names to use.
    use_raw : bool, optional
        Whether to use `adata.raw` values. Defaults to True.

    Returns
    -------
    anndata.AnnData
        Updated AnnData object with sex assignment columns in `.obs`.

    Examples
    --------
    >>> import scanpy as sc
    >>> import scanpy_plus as scp
    >>> scp.pp.assign_sex(adata, ygenes=['Ddx3y', 'Kdm5d'], use_raw=True)
    """

    if use_raw == True:
        xist_expr = data.raw[:,'Xist'].X.copy()
        if isinstance(xist_expr,scipy.sparse.csr.csr_matrix):
            xist_expr = xist_expr.toarray()
        data.obs['xist_logn'] = xist_expr[:,0]

        ygenes = data.raw.var.index.isin(ygenes)
        yexpr = data.raw.X[:,ygenes]
        if isinstance(yexpr, scipy.sparse.csr.csr_matrix):
            yexpr = yexpr.toarray()
    else:
        xist_logn = data[:,'Xist'].X.copy()
        if isinstance(xist_logn, scipy.sparse.csr.csr_matrix):
            xist_logn = xist_logn.toarray()
        data.obs['xist_logn'] = xist_logn[:,0]

        ygenes = data.var.index.isin(ygenes)
        yexpr = data.X[:,ygenes]
        if isinstance(yexpr, scipy.sparse.csr.csr_matrix):
            yexpr = yexpr.toarray()

    yexpr = yexpr.mean(axis = 1)
    data.obs['Ygene_logn'] = yexpr

    #plt.scatter(data.obs.xist_logn, data.obs.Ygene_logn, alpha = 0.05)
    #plt.show()

    data.obs['xist_bin'] = data.obs['xist_logn'] > 0
    data.obs['Ygene_bin'] = data.obs['Ygene_logn'] > 0
    temp = pd.crosstab(index = data.obs.xist_bin, columns = data.obs.Ygene_bin)
    sns.heatmap(temp, annot = True)
    print(temp)
    
    data.obs['sex'] = 'UKN'
    data.obs.loc[data.obs.xist_bin & ~data.obs.Ygene_bin , data.obs.columns == 'sex'] = 'F'
    data.obs.loc[data.obs.xist_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'MF'
    data.obs.loc[~data.obs.xist_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'M'
    
