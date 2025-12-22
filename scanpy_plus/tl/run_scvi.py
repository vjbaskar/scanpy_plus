import anndata as ad 
from typing import Union,List
import scanpy as sc
import matplotlib.pyplot as plt
from loguru import logger

def download_cc():
    """
    Download cellcycle genes
    """
    import requests
    import os

    # The official URL for the gene list
    url = "https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"

    # The desired local filename
    filename = "regev_lab_cell_cycle_genes.txt"

    # Check if the file already exists before downloading
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for bad status codes
            with open(filename, 'w') as f:
                f.write(response.text)
                print(f"Successfully downloaded {filename}.")
        except requests.exceptions.RequestException as e:
            print(f"Error downloading the file: {e}")
    else:
        print(f"{filename} already exists. Skipping download.")


def download_cc():
    """
    Download cellcycle genes
    """
    import requests
    import os

    # The official URL for the gene list
    url = "https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"

    # The desired local filename
    filename = "regev_lab_cell_cycle_genes.txt"

    # Check if the file already exists before downloading
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for bad status codes
            with open(filename, 'w') as f:
                f.write(response.text)
            print(f"Successfully downloaded {filename}.")
        except requests.exceptions.RequestException as e:
            print(f"Error downloading the file: {e}")
    else:
        print(f"{filename} already exists. Skipping download.")

def clean_genes(data, additional_genes_to_exclude =[]):
    """
    Clean up gene list. Remove cc, hb, mt, ribo and custom list of genes
    """
    from loguru import logger 
    import pandas as pd
    download_cc()
    cc_genes_csv=pd.read_csv('regev_lab_cell_cycle_genes.txt',  names=["gene_ids"], skiprows=1)
    cc_genes_csv = cc_genes_csv["gene_ids"]
    cc_genes_csv = list(cc_genes_csv)

    # Mark MT/ribo/Hb/cell cycle genes
    data.var['mt'] = data.var_names.str.startswith('MT-')  
    data.var["ribo"] = data.var_names.str.startswith(("RPS", "RPL"))
    data.var["hb"] = data.var_names.str.contains(("^HB[^(P)]")) 
    #data.var["hb"] = data.var_names.str.startswith(("HBA1", "HBA2", "HBB", "HBD","HBM", "HBZ", "HBG1", "HBG2", "HBQ1"))
    data.var["cc"] = data.var_names.isin(cc_genes_csv)
    mask_to_exclude = (
        data.var.cc | 
        data.var.hb | 
        data.var.mt |
        data.var.ribo |
        data.var.index.isin(additional_genes_to_exclude)
    )
    logger.info(f"Genes before removal: {data.shape[1]}")
    mask_to_include = ~mask_to_exclude
    data  = data[:, mask_to_include].copy()
    logger.info(f"Genes after removal: {data.shape[1]}")
    return data

def scvi_plot(model):
    """
    Plots to help if scvi has learnt
    If training and validation curves converge → good fit

    If validation loss diverges → possible overfitting or learning-rate issue

    If loss is noisy or doesn’t decrease → model may not be converging (try smaller LR, more epochs, or normalize input better)
    """
    import matplotlib.pyplot as plt

    plt.plot(model.history["elbo_train"], label="train")
    plt.plot(model.history["elbo_validation"], label="validation")
    plt.xlabel("Epochs")
    plt.ylabel("ELBO")
    plt.legend()
    plt.show()


def run_scvi(adata_hvg, 
    batch_key,
    raw_count_layer = 'counts',
    clean_genes_before_integration = True,
    n_latent=10, 
    n_layers=1,
    max_epochs=10,
    batch_size=512,
    categorical_covariate_keys=[], 
    continuous_covariate_keys=[],
    dispersion = 'gene-batch',

    latent_key='X_scvi',
    **kwargs
    ):
    """
    Run scvi
    """
    from loguru import logger 
    import scvi
    scvi.settings.seed = 42
    if clean_genes_before_integration:
        adata_hvg = clean_genes(adata_hvg)
    
    logger.info("Init model")
    scvi.model.SCVI.setup_anndata(adata_hvg, 
                                      layer=raw_count_layer,
                                      categorical_covariate_keys=categorical_covariate_keys,
                                      continuous_covariate_keys = continuous_covariate_keys,
                                      batch_key=batch_key,
                                      #                                labels_key="broad_annotation",
                                      # unlabeled_category="New/unlabelled/excluded"
                                     )
    model = scvi.model.SCVI(adata_hvg, dispersion= dispersion, n_latent = n_latent, n_layers = n_layers)
    logger.info("Started training model")

    #train_kwargs = {k: v for k, v in kwargs.items() if k in vae.train.__code__.co_varnames + run_scvi.train.Trainer.__init__.__code__.co_varnames}
    #vae.train(**train_kwargs)
    logger.info("Training model")
    model.train(max_epochs=max_epochs,             
                early_stopping=True,
                # accelerator='gpu',
                early_stopping_patience=5, #use_gpu =True, 
                batch_size=batch_size)
    logger.info("Plotting ELBO loss")
    scvi_plot(model)
    logger.info(f"Getting obsm updated. Your scvi model will be stored in {latent_key}")
    adata_hvg.obsm[latent_key] = model.get_latent_representation()
    return adata_hvg, model

