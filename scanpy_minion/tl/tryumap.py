from ..globimport import *
import time


def tryumap(adata, reqCols = ['total_counts'], regress_vars = ['S_score','G2M_score'], 
            regress = 0,
            scale = True,
            n_jobs = 4, 
            use_highly_variable = True, 
            ncols = 4,
            n_pcs = 50,
            n_neighbors = 15, 
            knn= True, 
            max_value_scale = 8, 
            save_prefix = None, 
            reset_fig_params = True, **kwargs):
    """
    Generates UMAP from anndata
    ---------------------------
    
    Input MUST be LOGNORM.
    If adata.X is scaled then use scale = False option
    
    
    Try generating umap, plot reqCols.
    data: anndata. Expects anndata with logcounts as X.
    regress: bool. If you want to regress out `regressed_vars` then set it to 1
    scale: True. Scale the data. This assumes that X is lognormed
    n_jobs: int. n procs for regress job
    use_highly_variable: bool. Used by sc.tl.pca
    ncols: int. Cols per row in UMAP. Used by sc.pl.umap
    n_pcs: int. Number of pcs to be used for sc.pp.neigbors
    n_neighbors, knn: int. Use by sc.pp.neighbours
    max_value_scale: float. Used by sc.pp.scale
    save_prefix: str. Write h5ad file as savedir_time.h5ad
    kwargs: UMAP arguments
    
    returns: anndata. 
    
    Example: tryumap(adata.copy(), reqCols = ['doublet_score', 'leiden'], regress = 0, n_jobs = 8, savedir="outputdir/prefixname") 
    No regression, scale data. Here n_jobs var will not be used as no sc.pp.regress_out will be run.
    
    """
    start = time.time()
    data = adata.copy(); print('Copy data')
    if reset_fig_params == True:
        sc.set_figure_params(scanpy=True)
    if regress == 1:
        print("Regressing")
        sc.pp.regress_out(data, keys = regress_vars, n_jobs= n_jobs)
    
    if use_highly_variable == True:
        try:
            print("HVG")
            print(data.var['highly_variable'].value_counts())
        except:
            print("Running sc.pp.highly_variable: highly_variable not found")
            sc.pp.highly_variable_genes(data) ; print("HVG")
            print(data.var['highly_variable'].value_counts())
    else:
        print("You have chosen to not run for highly variable")
        print("The PCA will be generated for vars, ", data.shape[1])
        
    if scale is True:
        data.raw = data
        sc.pp.scale(data, max_value = max_value_scale); print("Scaling")

    
    #try:
    sc.pp.pca(data, use_highly_variable=use_highly_variable) ; print('PCA')
    #except ValueError:
    #    print("Running sc.pp.highly_variable: highly_variable not found")
    #    sc.pp.highly_variable_genes(data) ; print("HVG")
    #    print(data.var['highly_variable'].value_counts())
    #    sc.pp.pca(data, use_highly_variable=use_highly_variable) ; print('PCA')
        
    sc.pl.pca_overview(data, color = reqCols)
    sc.pp.neighbors(data, n_neighbors = n_neighbors, knn = knn,  n_pcs = n_pcs  ) ; print("NN")
    sc.tl.umap(data) ; print("UMAP")
    #sc.tl.leiden(data) ; print("Leiden")
    sc.pl.umap(data, color = reqCols, ncols = ncols, **kwargs)
    if save_prefix is not None:
        from datetime import datetime
        now = datetime.now()
        nt = now.strftime("%Y_%m_%d_%H_%M_%S")
        f = "_".join([save_prefix, nt])
        f = f+".h5ad"
        print("Writing h5ad to ", f)
        data.write(f)
    print(f"Total time lapsed (sec) = {time.time() - start}")
    return data
