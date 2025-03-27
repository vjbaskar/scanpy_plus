def liana_cc_rankaggr(adata):

    """
    This function calculates consensus ligand-receptor predictions  of individual methods. This is done by ranking and aggregating (RRA) the ligand-receptor interaction predictions from all methods.
    adata: anndata obj
    
    return: anndata with liana ligand-receptor interactions consensus results
    """
	li.mt.rank_aggregate(adata_filtered_dr,use_raw=False,
                     groupby='celltype',
                     resource_name='consensus',
                     expr_prop=0.1
                     verbose=True, **kwargs)
	rank_aggregate.describe()

