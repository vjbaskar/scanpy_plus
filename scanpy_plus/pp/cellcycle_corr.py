from ..globimport import *
import time

def cellcycle_corr(data, cell_cycle_genes = None, sim_cutoff = 0.2):
    
    """
    Finds genes highly correlated with cell cycle genes.
    
    data: anndata. X MUST be LOGNORM
    cell_cycle_genes: list of cell cycle genes(1)
    sim_cutoff: Similarity cutoff
    
    returns: list(list) gene_correlation = matrix
                        genes_to_retain = genes to retain
    (1) cell_cycle_genes = [x.strip().capitalize() for x in open('public_data/regev_lab_cell_cycle_genes.txt')]
        For human. no need to capitalize()
    
    """
    #print(cell_cycle_genes)
    if cell_cycle_genes == None:
        print("No cell cycle gene")
        return [0,0]
    start = time.time()
    temp = sc.pp.scale(data, max_value=10, copy=True)
    gene_correlation = np.corrcoef(temp.X.transpose())
    gene_correlation_cell_cycle = gene_correlation[[i for i, gene in enumerate(temp.var_names) if gene in cell_cycle_genes], :]
    high_correlation_gene_filter = np.amax(gene_correlation_cell_cycle, axis=0) < sim_cutoff
    cell_cycle_removed_genes = temp.var_names[np.invert(high_correlation_gene_filter)]
    genes_to_retain = temp.var_names[high_correlation_gene_filter]
    #temp = temp[:, high_correlation_gene_filter]
    #return temp
    print(f"Total time taken (sec) = {time.time() - start}")
    return gene_correlation, genes_to_retain
    
    
