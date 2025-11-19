from typing import Optional
import anndata as an

def cellbender(
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
    import logging
    import sys
    import anndata as an
    import scanpy as sc
    from typing import Optional
    logging.basicConfig(level=logging.INFO)

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

