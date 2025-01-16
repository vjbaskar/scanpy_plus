from ..globimport import *

markers = dict(
    HSC=['Procr', 'Mllt3', 'Hoxb5', 'Mecom', 'Ly6a'],
    Ery=['Klf1', 'Gata1', 'Hba-a1', 'Epor'],
    MC=['Cma1', 'Gzmb', 'Kit'],
    Bas=['Prss34', 'Mcpt8'],
    Eos=['Prg2', 'Prg3'],
    Mk=['Pf4', 'Ppbp', 'Vwf'],
    Mo_DC=['Csf1r', 'Ctss', 'Ms4a6c', 'Irf8','Ccr2', 'F13a1', 'Mpo'],
    Neu=['Elane', 'Prtn3', 'Cebpe', 'Gfi1'],
    Ly=['Dntt', 'Il7r', 'Jchain', 'Ighm', 'Flt3', 'Rag2'],
    Bcell=['Cd79a', 'Cd79b', 'Cd19'],
    Tcell=['Cd3e', 'Bcl11b', 'Lck', 'Cd4', 'Cd8a'],
    InnateLymCells =['Gata3', 'Id2'],
    pDC=['Irf8', 'Siglech', 'Ly6d'],
    DC=['H2-Aa', 'Xcr1', 'Itgax'],
    Ifn_act=['Ifitm3', 'Irf7','Isg15'],
    MyoC1=['C1qa', 'C1qb', 'Mpo'],
    Nk=['Nkg7'],
    
              )

def haem_markers(adata, **kwargs):
    for k,v in markers.items():
        olap = list(set(adata.var_names) & set(v))
        if len(olap) > 0:
            sc.pl.umap(adata, color=olap, ncols=len(olap), title = list(map(lambda x: "/".join([k,x]), olap)), **kwargs)






