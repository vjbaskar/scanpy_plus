

def clean_batches(adata, batch_key, new_batch_key, min_cells=1000):
    """
    Clean up batches
    adata: anndata.AnnData
    batch_key: str
    new_batch_key: str
    min_cells: int  default 1000
    Returns:
    adata: anndata.AnnData
    new_batch_key: str
    Example:
    adata, new_batch_key = clean_baches(adata, batch_key, new_batch_key, min_cells=1000)
    """
    
    from loguru import logger
    import scanpy as sc
    import numpy as np
    # Creates a new covariate batch key
    adata.obs[new_batch_key] = adata.obs[batch_key].copy()

    batch_counts = adata.obs[new_batch_key].value_counts()
    small_batches = batch_counts[batch_counts < min_cells].index
    adata.obs[new_batch_key] = adata.obs[new_batch_key].astype(str)  # ensure string
    adata.obs[new_batch_key] = adata.obs[new_batch_key].replace(small_batches, "other")

    print("Number of batches before merging:", adata.obs[batch_key].nunique())

    print("Number of batches after merging:", adata.obs[new_batch_key].nunique())
    print(adata.obs[new_batch_key].value_counts().head())
    return adata, new_batch_key



def select_hvg(
    adata,
    ngenes,
    nbatches_min,
    batch_key
):
    """
    Select highly variable genes based on the number of batches and the number of genes
    adata.X must be log normalized
    ngenes: number of genes to select
    nbatches_min: minimum number of batches to select genes
    batch_key: batch key to use for selecting highly variable genes
    Returns:
    adata: anndata.AnnData
    nbatches_min: int minimum number of batches to select genes
    hvg_number: int number of highly variable genes
    label_dict: dict dictionary of labels for the number of batches
    label_dict2: dict dictionary of labels for the number of genes
    Example:
    hvg_var = select_hvg(adata, ngenes, nbatches_min, batch_key)
    """
    from loguru import logger
    import numpy as np
    import scanpy as sc
    
    logger.info("Cleaning batches: relegate small batches out")
    new_batch_key = batch_key + "_1"
    adata, new_batch_key = clean_batches(adata, batch_key, new_batch_key, min_cells=1000)
    
    logger.info(f"Selecting HVGs with {ngenes} genes and {nbatches_min} batches")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=ngenes,
        subset=False,
        batch_key=new_batch_key,
        check_values=False,
    )

    var_genes_all = adata.var.highly_variable
    var_genes_batch = adata.var.highly_variable_nbatches > nbatches_min
    var_select = adata.var.highly_variable_nbatches >= nbatches_min
    var_genes = var_select.index[var_select]

    hvg_number = len(var_genes)

    label_dict = adata.var['highly_variable_nbatches'].to_dict()
    label_dict2 = adata.var['highly_variable'].to_dict()


    best_nbatches_min = None
    closest_hvg_number = None
    closest_difference = float('inf')
    
    max_nbatches = len(adata.obs[batch_key].unique())
    logger.info(f"Max nbatches {max_nbatches}")    
    for nbatches_min in list(np.arange(1, max_nbatches, 1)):
        var_genes_batch = adata.var.highly_variable_nbatches > nbatches_min
        var_select = adata.var.highly_variable_nbatches >= nbatches_min
        var_genes = var_select.index[var_select]
        hvg_number = len(var_genes)

        difference = abs(hvg_number - ngenes)

        if difference < closest_difference:
            closest_difference = difference
            closest_hvg_number = hvg_number
            best_nbatches_min = nbatches_min

    nbatches_min = best_nbatches_min
    hvg_number = closest_hvg_number
    logger.info(f"selected {hvg_number} HVGs!")
    var_select = adata.var.highly_variable_nbatches >= nbatches_min
    logger.info(f"{hvg_number} selected -> {adata.shape}")
    return adata.var[var_select]

