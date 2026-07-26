Intermediate tutorial: plotting and data export
===============================================

Goal
----

Generate grouped plots and export/reload AnnData for reproducible sharing.

Inputs
------

- ``AnnData`` with categorical obs columns (for example ``condition`` and ``celltype``).

Steps
-----

.. code-block:: python

   import scanpy_plus as sp

   sp.pl.splitplot(adata, split_by="condition", colour_by="celltype")
   sp.pl.ridgeplot(adata, keys="log1p_total_counts", split_by="condition")

   sp.io.write_adata_split(adata, "analysis_export")
   reloaded = sp.io.read_adata_split("analysis_export", read_X=True)

Expected output and interpretation
----------------------------------

- Split and ridge plots enable side-by-side condition comparisons.
- ``analysis_export`` contains partitioned files suitable for sharing large objects.
- ``reloaded`` should match key annotations and structure from the source ``adata``.
