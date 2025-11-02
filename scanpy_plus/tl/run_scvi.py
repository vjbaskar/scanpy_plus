import anndata as ad 
from typing import Union,List
import scanpy as sc
import matplotlib.pyplot as plt


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


def run_scvi(adata_hvg,
             BATCH_KEY, 
             clean_genes = True,
             N_LATENT=10, 
             N_LAYERS=1,
             MAX_EPOCHS=10,
             BATCH_SIZE=512,
             CATEGORICAL_COV=[], 
             CONTINUOUS_COV=[],
             DISPERSION = 'gene-batch',
             **kwargs
            ):
    import scvi
    if clean_genes:
        adata_hvg = clean_genes(adata_hvg)
        scvi.model.SCVI.setup_anndata(adata_hvg, 
                                      layer="counts",
                                      categorical_covariate_keys=CATEGORICAL_COV,
                                      continuous_covariate_keys = CONTINUOUS_COV,
                                      batch_key=BATCH_KEY,
                                      #                                labels_key="broad_annotation",
                                      # unlabeled_category="New/unlabelled/excluded"
                                     )
        model = scvi.model.SCVI(adata_hvg, 
                                dispersion=DISPERSION,
                                n_latent = N_LATENT, 
                                n_layers = N_LAYERS,
                               )

    #train_kwargs = {k: v for k, v in kwargs.items() if k in vae.train.__code__.co_varnames + run_scvi.train.Trainer.__init__.__code__.co_varnames}
    #vae.train(**train_kwargs)
    model.train(max_epochs=MAX_EPOCHS,             
                early_stopping=True,
                # accelerator='gpu',
                early_stopping_patience=5, #use_gpu =True, 
                batch_size=BATCH_SIZE)


    print("model trained")
    return adata_hvg, model
