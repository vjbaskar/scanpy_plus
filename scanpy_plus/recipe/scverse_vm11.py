
import anndata as an
import numpy as np
import scanpy as sc
import sys
import logging
from typing import Optional


def is_outlier(adata, metric: str, nmads: int):
    from scipy.stats import median_abs_deviation
    import numpy as np

    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier



def run_cellbender_remove_background(
    adata: an.AnnData,
    #output_adata: an.AnnData,
    expected_cells: Optional[int] = None,
    total_droplets_included: Optional[int] = None,
    epochs: Optional[int] = None,
    use_cuda: bool = False,
    fpr: Optional[float] = None,
    extra_args: list = None
):
    """
    Run CellBender remove-background via base_cli.

    Parameters:
        input_path (str): Path to raw feature-barcode matrix HDF5 file.
        output_path (str): Output path for cleaned file.
        expected_cells (int): Expected number of cells.
        total_droplets_included (int): Number of droplets to include.
        epochs (int): Number of training epochs.
        use_cuda (bool): Whether to use CUDA (GPU).
        fpr (float): False positive rate for cell calling.
        extra_args (list): Any additional command-line flags (as a list).
    """
    from cellbender import base_cli
    import tempfile

    cellbender = tempfile.mkdtemp()
    cellbender_input = f"{cellbender}/cellbender_input.h5ad"
    cellbender_output = f"{cellbender}/cellbender_output.h5ad"
    adata.write_h5ad(cellbender_input)
    args = [
        "remove-background",
        "--input", cellbender_input,
        "--output", cellbender_output
    ]
    if expected_cells is not None:
        args.append("--expected-cells")
        args.append(str(expected_cells))
    if total_droplets_included is not None:
        args.append("--total-droplets-included")
        args.append(str(total_droplets_included))
    if epochs is not None:
        args.append("--epochs")
        args.append(str(epochs))
    if fpr is not None:
        args.append("--fpr")
        args.append(str(fpr))
    if not use_cuda:
        pass
    else:
        args.append("--cuda")

    if extra_args:
        args += extra_args

    logging.info("Running CellBender with arguments: %s", args)

    # Fake command-line input
    sys.argv = ["cellbender"] + args


    try:
        base_cli.main()
        
    except SystemExit as e:
        if e.code != 0:
            logging.error("CellBender failed with exit code: %s", e.code)
            raise RuntimeError(f"CellBender exited with code {e.code}")
        else:
            logging.info("CellBender completed successfully.")
            output_adata = an.read_h5ad(cellbender_output)
            return output_adata


def scverse_vm11(
    adata: an.AnnData,
    nmads: float = 5,
    mt_nmads: float = 3,
    mt_threshold: float = 10,
    cellbender: bool = False,
    hvg_flavour: str = "seurat_v3",
    n_top_genes: int = 2000,
    n_comps: int = 50,
    n_neighbours: int = 10,
    min_dist: float = 0.5,
    spread: float = 1.0,
    cluster_resolution: float = 0.5,
    scale = False
):
    """
    Perform standard single-cell RNA-seq preprocessing and clustering using the scverse stack.

    This pipeline includes filtering, highly variable gene selection, dimensionality reduction,
    neighborhood graph construction, UMAP embedding, and Leiden clustering.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.
    nmads : float, default 5
        Threshold (in MADs) for general outlier filtering (e.g. low quality cells).
    mt_nmads : float, default 3
        Threshold (in MADs) for filtering based on mitochondrial gene expression.
    mt_threshold : float, default 10
        Absolute threshold (% of mitochondrial gene expression) beyond which cells are filtered.
    cellbender : bool, default False
        Whether to skip filtering steps assuming CellBender has already removed ambient RNA.
    hvg_flavour : str, default "seurat_v3"
        Method for selecting highly variable genes. Common options: "seurat_v3", "cell_ranger".
    n_top_genes : int, default 2000
        Number of highly variable genes to retain.
    n_comps : int, default 50
        Number of principal components to compute.
    n_neighbours : int, default 10
        Number of neighbors to compute in the neighborhood graph.
    min_dist : float, default 0.5
        Minimum distance parameter for UMAP.
    spread : float, default 1.0
        Spread parameter for UMAP.
    cluster_resolution : float, default 0.5
        Resolution parameter for Leiden clustering.
    scale : bool, default False
        Whether to scale the data to unit variance before PCA.

    Returns
    -------
    None
        The function modifies the `adata` object in-place with the results of preprocessing,
        dimensionality reduction, and clustering (e.g., adding `.obs['leiden']`, `.obsm['X_umap']`, etc.).
    
    Example
    -------
    >>> import scanpy as sc
    >>> from scanpy_plus.recipe import scverse_vm11
    >>> adata = sc.read("path_to_your_data.h5ad")
    >>> adata = scverse_vm11(adata, nmads=5, mt_nmads=3, mt_threshold=10, cellbender=False,
    ...                     hvg_flavour="seurat_v3", n_top_genes=2000, n_comps=50,
    ...                     n_neighbours=10, min_dist=0.5, spread=1.0, cluster_resolution=0.5, scale=False)
    """

    import numpy as np
    import scanpy as sc
    import seaborn as sns
    from scipy.stats import median_abs_deviation
    import tempfile
    from cellbender.remove_background import run

    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=300,
        facecolor="white",
        frameon=False,
    )

    if cellbender:
        # Run CellBender remove-background
        if adata.n_obs < 20000:
            logging.warning("Cellbended needs raw data. Looks like there is less cells. Skipping cellbender")
        else:
            adata = run_cellbender_remove_background(adata)

    # 1. Preprocess the data
    adata.raw = adata 
    adata.layers["counts"] = adata.X.copy()
    ## Compute genes
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    ## Compute QC metrics
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

    ## Plot QC metrics
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    # sc.pl.violin(adata, 'total_counts')
    p2 = sc.pl.violin(adata, "pct_counts_mt")
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", nmads)
    | is_outlier(adata, "log1p_n_genes_by_counts", nmads)
    | is_outlier(adata, "pct_counts_in_top_20_genes", nmads)
    )
    print(adata.obs.outlier.value_counts())

    # mt filter
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_nmads) | (adata.obs["pct_counts_mt"] > mt_threshold)
    adata.obs.mt_outlier.value_counts()

    # Filter out low quality cells
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

    # 2. Flag doublets
    sc.pp.scrublet(adata)

    # 3. Normalize the data
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    del scales_counts

    # 4. select highly variable genes
    # Change the adata options for different flavours
    if hvg_flavour == "seurat_v3":
       sc.pp.highly_variable_genes(adata,
                                n_top_genes=n_top_genes,
                                flavor=hvg_flavour, subset=True)
       
    else:
        raise ValueError("Unknown flavour")


    ax = sns.scatterplot(
        data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
    )
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    plt.show()

    if scale:
        sc.pp.scale(adata, max_value=10, zero_center=True)
    
    # 5. PCA
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)

    # 6. Neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbours, n_pcs=n_comps)

    # 7. UMAP
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=42)

    # 8. Clustering
    sc.tl.leiden(adata, resolution=cluster_resolution, random_state=42)
    return adata
