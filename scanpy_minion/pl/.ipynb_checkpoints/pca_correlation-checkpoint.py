from ..globimport import *

def pca_correlation(data, obsm_pca = 'X_pca', fill_na = 0, figsize =(25,25), fontsize = 0.75, clustermap = True, **kwargs):
    
    """
    Calculates correlation of each PCs to obs data.
    Automatically removes columns that are string or category.
    
    Parameters:
    data: anndata object
    obsm_pca: The name of obsm_pca to consider. eg. X_pca, X_umap, X_pca_wnn
    
    Returns: pandas dataframe
    
    """
    
    # Generate obs + pca object
    print('Concatenating obs and pca')
    obs_pca = data.obs.copy()
    lsize = data.obsm[obsm_pca].shape[1]
    print(lsize)
    temp = pd.DataFrame(data.obsm[obsm_pca], index=data.obs_names, columns=['PCA'+str(i) for i in range(lsize)])
    obs_pca = pd.merge(obs_pca, temp, how='left', left_index=True, right_index=True)
     
    # Retain only columns that are float or int
    print("Removing columns with str or cat dtypes")
    obs_pca = obs_pca.select_dtypes(include=['float', 'int', 'float32', 'float64', 'int32', 'int64'])

    # Calculating correlation matrix
    print("Computing correlation matrix")
    cols = obs_pca.columns
    corr_matx = np.zeros(shape = (len(cols), len(cols)))
    a=0
    b=0
    for i in obs_pca.columns:
        for j in obs_pca.columns:
            x = np.array(obs_pca[i])
            y = np.array(obs_pca[j])
            corr_matx[a,b] = np.corrcoef(x,y)[0,1]
            b = b + 1
        a = a + 1
        b = 0
        print(".", end="")
    print()
    
    print("Plotting corr matrix")
    corr_df = pd.DataFrame(corr_matx, index=cols, columns=cols)
    corr_df[corr_df.isna()] = fill_na
    plt.rcParams['figure.figsize'] = figsize
    sns.set(font_scale=fontsize)
    if clustermap == True:
        sns.clustermap(corr_df, cmap= 'vlag', yticklabels=True, xticklabels=True, **kwargs)
    print("sns.clustermap(corr_df, cmap= 'vlag', yticklabels=True, xticklabels=True)")
    sc.set_figure_params(scanpy=True, color_map='GnBu')
    return corr_df