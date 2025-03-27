def liana_cc_dotplot(adata, 
                     source_celltype:str list, 
                     target_celltype: str list, **kwargs):
        
    """
    Plot LIANA results
    adata: anndata obj
    source_celltype: celltype(s) to be used as ligands (string)
    target_celltype: celltype(s) to be used as ligands (string)
    
    return: dotplot with LR results
    """


    li.pl.dotplot(adata = adata_filtered_dr,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=source_celltype,
              target_labels=target_celltype,
              top_n=30,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(8, 9)
             )
