from ..globimport import *
import scipy
import seaborn as sns
def assign_sex(data, ygenes, use_raw = True):
    """
    Assigns sex based on Xist and Y-chromosome gene expression. 
    Expects log-normalised input.
    Thanks to Iwo for the code.
    Adds these cols to obs: xist_logn, Ygene_logn, xist_bin, Ygene_bin, sex
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
    
