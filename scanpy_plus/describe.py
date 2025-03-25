from .globimport import *
import numpy as np
def describe(data):
    """
    Descibes the X of anndata. Useful for checking if X is normalised, raw etc
    data: anndata obj
    returns: None
    """
    X = data.X
    dmin = X.min()
    dmax = X.max()
    cellsum = np.sum(data.X[1:100,], axis = 1)
    cellmin = cellsum.min()
    cellmax = cellsum.max()
    cellmean = cellsum.mean()
    celldiff = cellmax - cellmin
    ncells, ngenes = data.shape
    
    #print("****", end = " ")
    print("****", end = " ")
    if dmin == 0:
        if dmax > 1000:
            if (celldiff > 1.4* cellmean):
                print("Highly likely data is raw")
            else:
                print("Highly likely data is normalised")
        if dmax < 50:
            print("Highly likely data is logged. Probably lognorm-ed")
    else:
        print("Highly likely data is scaled or regressed")
    print("****")
    print("-----")
    print(f"(#) Cells = {ncells}", end = " | ")
    print(f"Genes = {ngenes}", end = "\n")
    print(f"(X) min = {dmin}", end = " | ")
    print(f"max = {dmax}", end = " \n")
    print("-----")
    print("(head 100 cells) umi_sum min = ", cellmin, end = " | ")
    print("umi_sum max = ", cellmax, end = " | ")
    print("cellmean =  ", cellmean, end = "\n")
    print("-----")
    print("Obs: ", data.obs.columns.tolist())
    print("Var: ", data.var.columns.tolist())

def describe(data):
    """
    Descibes the X of anndata. Useful for checking if X is normalised, raw etc
    data: anndata obj
    returns: None
    """
    X = data.X
    dmin = X.min()
    dmax = X.max()
    cellsum = np.sum(data.X[1:100,], axis = 1)
    cellmin = cellsum.min()
    cellmax = cellsum.max()
    cellmean = cellsum.mean()
    celldiff = cellmax - cellmin
    ncells, ngenes = data.shape
    
    print("****", end = " ")
    if dmin == 0:
        if dmax > 1000:
            if (celldiff > 1.4* cellmean):
                print("Highly likely data is raw")
            else:
                print("Highly likely data is normalised")
        if dmax < 50:
            print("Highly likely data is logged. Probably lognorm-ed")
    else:
        print("Highly likely data is scaled or regressed")
    print("-----")
    print(f"(X) min = {dmin}", end = " | ")
    print(f"max = {dmax}", end = " \n")
    print("-----")
    print("(head 100 cells) umi_sum min = ", cellmin, end = " | ")
    print("umi_sum max = ", cellmax, end = " | ")
    print("cellmean =  ", cellmean, end = "\n")
    print("-----")
    print(data)

anndata.AnnData.describe = describe
