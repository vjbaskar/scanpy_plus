Usage
=====

Assumptions
-----------

Examples below assume:

- You have an ``AnnData`` object called ``adata``.
- ``adata.raw`` is populated when using raw-expression helpers.
- Common QC columns (for example ``total_counts`` and ``n_genes_by_counts``) are available when referenced.

Basic workflow
--------------

.. code-block:: python

   import scanpy as sc
   import scanpy_plus as sp

   adata = sc.read_h5ad("example.h5ad")

   # Optional: restore matrix from adata.raw
   adata = sp.get_from_raw(adata)

   # QC helpers
   adata.obs["is_outlier"] = sp.is_outlier(adata, metric="total_counts", nmads=5)

   # Preprocessing helper
   adata = sp.pp.assign_sex(adata)

   # Dimension reduction helper
   adata = sp.tl.tryumap(adata)

   # Interactive UMAP export
   sp.pl.plotly_umap(
       adata,
       color_col="leiden",
       obs_cols=["sample", "condition"],
       write_to_file="umap.html",
   )

Common workflows
----------------

1. Batch-aware HVG selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   adata = sp.pp.clean_batches(adata, batch_key="batch", new_batch_key="batch_clean")
   adata_hvg = sp.pp.select_hvg(adata, batch_key="batch_clean")

Expected output: ``adata_hvg.var`` contains the HVG mask for downstream modeling.

2. SCVI-ready analysis
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   adata_hvg = sp.tl.clean_genes(adata)
   model = sp.tl.run_scvi(adata_hvg)

Expected output: a trained SCVI model and latent embedding stored in ``adata_hvg.obsm``.

3. Split and restore AnnData objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   sp.io.write_adata_split(adata, "adata_parts")
   adata_reloaded = sp.io.read_adata_split("adata_parts", read_X=True)

Expected output: a reconstructed ``AnnData`` object that can be re-used in later sessions.
