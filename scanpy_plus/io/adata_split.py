## Saving anndata

def write_adata_split(adata, output_folder):
    """
    Splits and saves anndata
    adata: annData [req]
    output_folder: str. Folder to save [ref]
    returns: None
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from scipy import sparse
    import os
    
    #output_folder = "adata_split"

    outdir = Path(output_folder)
    (outdir / "layers").mkdir(parents=True, exist_ok=True)
  
    X = adata.X
    cells = adata.obs_names
    genes = adata.var_names

    print(f"Writing X to {outdir}")
    if sparse.issparse(X):
        sparse.save_npz(outdir / "X.npz", X)
    else:
        np.save(outdir / "X.npy", X)
        
    print(f"Cells and genes")

    cells.to_series().to_csv(outdir / "cell_names.csv", index=False, header=["cell_names"])
    genes.to_series().to_csv(outdir / "var_names.csv", index=False, header=["var_names"])
    
    print(f"Saved X ({'sparse' if sparse.issparse(X) else 'dense'}) and names :: {outdir}")
    
    adata.obs.to_csv(outdir / "obs.csv")
    adata.var.to_csv(outdir / "var.csv")
    
    for layer_name, layer_data in adata.layers.items():
        layer_dir = outdir / "layers" # / layer_name
        #layer_dir.mkdir(parents=True, exist_ok=True)
    
        if sparse.issparse(layer_data):
            sparse.save_npz(layer_dir / f"{layer_name}.npz", layer_data)
        else:
            np.save(layer_dir / f"{layer_name}.npy", layer_data)
    
        print(f"Saved layer '{layer_name}' :: {layer_dir}")
    

    # Save anndata without X
    del adata.X
    adata.write_h5ad(outdir / "adata_metadata.h5ad", compression="gzip")
    print("Saved metadata-only AnnData :: adata_metadata.h5ad")

def read_adata_split(input_folder, read_X = False, read_layers = False):
    """
    Read anndata that is split.
    input_folder: str. Folder containing split anndata [req]
    read_X: bool. Whether to read the X matrix [ref]
    read_layers: bool. Whether to read the layers [ref]
    
    returns:
    anndata.AnnData - Reconstructed AnnData object
    """

    import anndata as ad
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from scipy import sparse
    import os
    from loguru import logger
    #input_folder = "adata_split"
    indir = Path(input_folder)

    logger.info(f"Reading AnnData from {indir}") 
    # Load metadata-only AnnData
    adata_meta = ad.read_h5ad(indir / "adata_metadata.h5ad")
    

        
    
    # Load cells and genes
    logger.info(f"Loading cell and gene names from {indir}")
    cells = pd.read_csv(indir / "cell_names.csv")["cell_names"]
    genes = pd.read_csv(indir / "var_names.csv")["var_names"]
    
    # Build AnnData
    adata_reconstructed = adata_meta
    
    # Load X
    if read_X:
        X_path = indir / "X.npz" if (indir / "X.npz").exists() else indir / "X.npy"
        logger.info(f"Loading X from {X_path}")
        try:
            X = sparse.load_npz(X_path)
        except:
            arr = np.load(X_path, allow_pickle=True)
            if isinstance(arr, np.lib.npyio.NpzFile):
                # Prefer key 'X', but fall back to first available
                key = "X" if "X" in arr.files else arr.files[0]
                X = arr[key]
            else:
                X = arr
        adata_reconstructed.X = X

    # Load layers
    if read_layers:
        
        layer_dir = indir / "layers"
        for layer_name in os.listdir(layer_dir):
            layer_name = layer_name.split(".")[0]  # remove extension
            layer_path = layer_dir  / f"{layer_name}.npz"
            logger.info(f"Loading layers from {layer_path}")
            if layer_path.exists():
                adata_reconstructed.layers[layer_name] = sparse.load_npz(layer_path)
            else:
                layer_path = layer_dir / f"{layer_name}.npy"
                adata_reconstructed.layers[layer_name] = np.load(layer_path)

        logger.info("Reconstructed full AnnData")
    return adata_reconstructed

