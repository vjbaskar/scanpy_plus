Intermediate tutorial: SCVI pipeline
====================================

Goal
----

Prepare data, run SCVI, and inspect downstream embeddings.

Inputs
------

- ``AnnData`` object with count-like data and batch metadata.

Steps
-----

.. code-block:: python

   import scanpy_plus as sp

   adata_hvg = sp.pp.select_hvg(adata, batch_key="batch")
   adata_hvg = sp.tl.clean_genes(adata_hvg)

   model = sp.tl.run_scvi(adata_hvg)

   # Visualize latent-space structure
   sp.pl.plotly_umap(
       adata_hvg,
       color_col="leiden",
       obs_cols=["batch"],
       write_to_file="scvi_umap.html",
   )

Expected output and interpretation
----------------------------------

- ``model`` is a fitted SCVI model for the selected features.
- The output HTML UMAP allows interactive inspection of batch mixing and cluster structure.
