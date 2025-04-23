import anndata as ad 
from typing import Union,List

# Define function for running scVI (from Veronika)
def run_scvi(adata : ad.AnnData, 
             layer_raw : str = 'counts',
             include_genes : list = [], 
             exclude_cc_genes : bool = True, 
             exclude_mt_genes : bool =True, 
             exclude_vdjgenes : bool = True, 
             remove_cite : bool =False,
             
             batch_hv : str = "age_group", 
             hvg : int = 3000,
             batch_scvi: str = "sample",
             cat_cov_scvi: list = ["DonorID", "10X_version", "Sex", "Age_group"], 
             cont_cov_scvi: list = ["percent_mito"], 
             leiden_clustering : float = 1.0, 
             col_cell_type : Union[str, List] = 'cell_type_level_4', 
             fig_prefix : str = 'scviEmbed',
             fig_dir : str = None,
             
             **kwargs):
    import scanpy as sc
    import scvi
    from typing import Union,List
    import matplotlib.pyplot as plt

   
    
    """Wrapper function for running scVI on anndata object.
    
    Parameters
    ----------
    adata : ad.AnnData object
        Annotated data matrix.
    layer_raw : str, optional
        Layer of adata to use for scVI, by default 'counts'.
    include_genes : list, optional
        List of genes to include in scVI, by default [].
    exclude_cc_genes : bool, optional
        Whether to exclude cell cycle genes, by default True.
    exclude_mt_genes : bool, optional
        Whether to exclude mitochondrial genes, by default True.
    exclude_vdjgenes : bool, optional
        Whether to exclude VDJ genes, by default True.
    remove_cite : bool, optional
        Whether to remove CITEseq genes, by default False.
    batch_hv : str, optional
        Batch key for highly variable gene selection, by default "age_group".
    hvg : int, optional
        Number of highly variable genes to select, by default 3000.
    batch_scvi : str, optional
        Batch key for scVI, by default "sample".
    cat_cov_scvi : list, optional
        Categorical covariates for scVI, by default ["DonorID", "10X_version", "Sex", "Age_group"]/
    cont_cov_scvi : list, optional
        Continuous covariates for scVI, by default ["percent_mito"].
    leiden_clustering : float, optional
        Resolution for leiden clustering, by default 1.0.
    col_cell_type : Union[str, List], optional
        Cell type annotation for UMAP, by default 'cell_type_level_4'.
    fig_dir : str | None, optional
        Directory to save UMAP plots, by default None.
    fig_prefix : str | None, optional
        Prefix for UMAP plots, by default 'scviEmbed'.
        
    Returns
    -------
    results : dict
        Dictionary containing the following
        - data : ad.AnnData object with scVI results
        - vae : scvi.model.SCVI object

    Example
    -------
    from scvi_wrapper import run_scvi

    object_version = f'v1_{today}'

    # Run scvi
    scvi_run = run_scvi(adata_concat,
                        layer_raw = 'X',
                        # Excluded genes
                        include_genes=[], exclude_cc_genes=True, exclude_mt_genes=True,
                        exclude_vdjgenes = True, remove_cite = False,
                        # Highly variable gene selection
                        batch_hv="chemistry", hvg = 5000, #span = 0.5,
                        # scVI
                        batch_scvi="sample_id",
                        cat_cov_scvi=["embryo", "chemistry", "sex", 'section'],
                        cont_cov_scvi=["percent_mito", 'n_genes'],
                        max_epochs=400, batch_size=2000, early_stopping = True, early_stopping_patience = 15, early_stopping_min_delta = 10.0,
                        plan_kwargs = {'lr': 0.001, 'reduce_lr_on_plateau' : True, 'lr_patience' : 10, 'lr_threshold' : 20},
                        n_layers = 3, n_latent = 30, dispersion = 'gene-batch',
                        # Leiden clustering
                        leiden_clustering = None, col_cell_type = ['cell_type_lvl1', 'cell_type_lvl5'],
                        fig_dir = f'{plots_path}/whole_embryo/preprocessing', fig_prefix = f'wholeEmbryo_allGEX_scvi_{object_version}')

    """
    
    # Copy adata
    counts = adata.layers[layer_raw] if layer_raw != 'X' else adata.X
    adata_scvi = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    
    # Remove excluded genes
    if remove_cite:
        print('removing CITEseq genes pre SCVI')
        adata_scvi = adata_scvi[:,~adata_scvi.var['cite']].copy() # remove cite genes
    gene_list = adata_scvi.var_names.tolist()
    if exclude_cc_genes:
        cell_cycle_genes = [x.strip() for x in open('/nfs/team205/vk8/processed_data/regev_lab_cell_cycle_genes.txt')]
        [gene_list.remove(i) for i in cell_cycle_genes if i in gene_list]
    if exclude_mt_genes:
        mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')]
        [gene_list.remove(i) for i in mt_genes if i in gene_list]
    if exclude_vdjgenes:
        import re
        [gene_list.remove(i) for i in gene_list if re.search('^TR[AB][VDJ]|^IG[HKL][VDJC]', i)]
    print('Removed excluded genes')
    
    # Select highly variable genes
    adata_scvi = adata_scvi[:,gene_list].copy()
    hvg_kwargs = {'adata' : adata_scvi, 'batch_key' : batch_hv, 'n_top_genes' : hvg, 'flavor' : 'seurat_v3'}
    hvg_kwargs.update({k:v for k,v in kwargs.items() if k in sc.pp.highly_variable_genes.__code__.co_varnames})
    sc.pp.highly_variable_genes(**hvg_kwargs)
    selected_genes = list(set(adata_scvi.var.loc[adata_scvi.var['highly_variable']].index.tolist() + include_genes))
    adata_scvi = adata_scvi[:, selected_genes].copy()
    print(f'Highly variable genes selected in total {adata_scvi.shape}')
    
    # Run scVI
    scvi.model.SCVI.setup_anndata(adata_scvi, batch_key=batch_scvi,
                                  categorical_covariate_keys=cat_cov_scvi,
                                  continuous_covariate_keys=cont_cov_scvi)
    scvi_kwargs = {k: v for k, v in kwargs.items() if k in scvi.model.SCVI.__init__.__code__.co_varnames + scvi.module.VAE.__init__.__code__.co_varnames}
    vae = scvi.model.SCVI(adata_scvi, **scvi_kwargs)
    train_kwargs = {k: v for k, v in kwargs.items() if k in vae.train.__code__.co_varnames + scvi.train.Trainer.__init__.__code__.co_varnames}
    vae.train(**train_kwargs)
    print('scVI model trained')
    
    # Perform UMAP
    adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
    sc.tl.umap(adata_scvi, min_dist = 0.1, spread = 1) 
    print("UMAP performed")
    
    # Perform clustering
    if bool(leiden_clustering) is True:
        if type(leiden_clustering) == float:
            leiden_clustering = [leiden_clustering]
        for res in leiden_clustering:
            sc.tl.leiden(adata_scvi, resolution = res, key_added = f"leiden_r{res}")
        print('Leiden clustering performed')
    
    # Construct new adata
    adata_raw_scvi = adata.copy()
    if bool(leiden_clustering) is True:
        if type(leiden_clustering) == float:
            leiden_clustering = [leiden_clustering]
        for res in leiden_clustering:
            adata_raw_scvi.obs[f"leiden_r{res}"] = adata_scvi.obs["leiden_r{res}"].copy()   
    adata_raw_scvi.obsm['X_scVI'] = adata_scvi.obsm['X_scVI'].copy()
    adata_raw_scvi.obsm['X_umap'] = adata_scvi.obsm['X_umap'].copy()
    adata_raw_scvi.obsp = adata_scvi.obsp.copy()
    adata_raw_scvi.uns = adata_scvi.uns.copy()
 
    # Plot UMAPs
    # Cell type annotations
    if type(col_cell_type) == str:
        col_cell_type = [col_cell_type] if col_cell_type in adata_raw_scvi.obs.columns else []
    elif type(col_cell_type) == list:
        col_cell_type = [c for c in col_cell_type if c in adata_raw_scvi.obs.columns]
    if any(col_cell_type):
            p = sc.pl.umap(
                adata_raw_scvi,
                color=col_cell_type, legend_loc = "on data", legend_fontsize = 4, frameon=False,
                ncols=2, return_fig = True, show = False)
            if fig_dir is not None:
                plt.savefig(f'{fig_dir}/{fig_prefix}_ctAnnotations_umap.png', dpi=300, bbox_inches='tight')
    # Covariates
    covs = []
    if cat_cov_scvi is not None:
        covs.append(cat_cov_scvi)
    if cont_cov_scvi is not None:
        covs.append(cont_cov_scvi)  
    if covs is not None:
        p = sc.pl.umap(
            adata_raw_scvi,
            color=cont_cov_scvi, frameon = False, ncols = 2,
            return_fig = True, show = False)
        if fig_dir is not None:
            plt.savefig(f'{fig_dir}/{fig_prefix}_covariates_umap.png', dpi=300, bbox_inches='tight')
   
    results = {}
    results['data'] = adata_raw_scvi
    results['vae'] = vae
    return (results)
