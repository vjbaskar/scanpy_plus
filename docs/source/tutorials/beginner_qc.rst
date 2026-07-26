Beginner tutorial: QC and first UMAP
====================================

Goal
----

Load an ``AnnData`` object, run basic QC-oriented helpers, and generate a first UMAP.

Inputs
------

- ``example.h5ad`` with an expression matrix and common obs QC columns.

Steps
-----

.. code-block:: python

   import scanpy as sc
   import scanpy_plus as sp

   adata = sc.read_h5ad("example.h5ad")
   adata.obs["is_outlier"] = sp.is_outlier(adata, metric="total_counts", nmads=5)

   adata = sp.pp.assign_sex(adata)
   adata = sp.tl.tryumap(adata)

   sp.pl.umap_allobs(adata, wanted_obs=["leiden", "sample"], cmap="viridis")

Expected output and interpretation
----------------------------------

- ``adata.obs['is_outlier']`` flags potential QC outliers.
- ``adata.obs`` gets sex-assignment annotations.
- UMAP visual output helps confirm whether major groups separate as expected.
