from ..globimport import *
import scipy
import seaborn as sns

def assign_sex(data, ygenes = None, use_raw = True, download_gencode = False):
    """
    Assigns sex based on XIST and Y-chromosome gene expression.

    This function expects log-normalised input. It adds the following columns to `adata.obs`:
    
    - `XIST_logn`: log-normalized expression of XIST
    - `Ygene_logn`: mean log-normalized expression of the specified Y-chromosome genes
    - `XIST_bin`: binary indicator for XIST expression
    - `Ygene_bin`: binary indicator for Y gene expression
    - `sex`: inferred sex ("male", "female", or "ambiguous")

    Parameters
    ----------
    data : anndata.AnnData
        Annotated data matrix.
    ygenes : list of str or None
        List of Y-chromosome gene names to use.
        If None: It uses default human chrY list from Ensembl.
    use_raw : bool, optional
        Whether to use `adata.raw` values. Defaults to True.
    download_gencode: bool
        Whether to get chrY genes from gencode and use it. Default False. [Time consuming]

    Returns
    -------
    anndata.AnnData
        Updated AnnData object with sex assignment columns in `.obs`.

    Examples
    --------
    >>> import scanpy as sc
    >>> import scanpy_plus as scp
    >>> scp.pp.assign_sex(adata, ygenes=['Ddx3y', 'Kdm5d'], use_raw=True)
    """
    import seaborn as sns
    import sys
    import scipy
    import pandas as pd
    from pathlib import Path


    if ygenes is None:
        if download_gencode:
            if not Path('gencode.basic.gz').is_file():
                print("Getting ygenes from gencode")
                subprocess.run(["wget", "-q",
                        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz",
                       "-O",
                       "gencode.basic.gz"], 
                     capture_output=True)

            full_command = """
            zcat gencode.basic.gz | \
            grep ^chrY | \
            awk ' $3 == "gene" ' | \
            cut -f 9 | \
            cut -d ";" -f 3 | \
            sed -e 's/ gene_name "//g' -e 's/"//g' | \
            sort -u
            """

            # Run the command using the shell
            # text=True decodes stdout/stderr as text
            # check=True will raise an exception if the command fails
            result = subprocess.run(full_command, shell=True, capture_output=True, text=True, check=True)
            ygenes = result.stdout.split("\n")
        else:
            print("Using default chrY gene set")
            ygenes = ["AKAP17A","AMELY","ASMT","ASMTL","BPY2",
                      "BPY2B","BPY2C","CD99","CDY1","CDY1B",
                      "CDY2A","CDY2B","CRLF2","CSF2RA","DAZ1","DAZ2",
                      "DAZ3","DAZ4","DDX3Y","DHRSX","EIF1AY","GTPBP6",
                      "HSFY1","HSFY2","IL3RA","IL9R","KDM5D","NLGN4Y",
                      "P2RY8","PCDH11Y","PLCXD1","PPP2R3B","RBMY1A1","RBMY1B",
                      "RBMY1D","RBMY1E","RBMY1F","RBMY1J","RPS4Y1","RPS4Y2","SHOX",
                      "SLC25A6","SRY","TBL1Y","TGIF2LY","TMSB4Y","TSPY1","TSPY10",
                      "TSPY2","TSPY3","TSPY4","TSPY8","TSPY9","USP9Y","UTY","VAMP7",
                      "VCY","VCY1B","WASH6P","ZBED1","ZFY"]

    if use_raw == True:
        XIST_expr = data.raw[:,'XIST'].X.copy()
        if isinstance(XIST_expr,scipy.sparse.csr.csr_matrix):
            XIST_expr = XIST_expr.toarray()
        data.obs['XIST_logn'] = XIST_expr[:,0]

        ygenes = data.raw.var.index.isin(ygenes)
        yexpr = data.raw.X[:,ygenes]
        if isinstance(yexpr, scipy.sparse.csr.csr_matrix):
            yexpr = yexpr.toarray()
    else:
        XIST_logn = data[:,'XIST'].X.copy()
        if isinstance(XIST_logn, scipy.sparse.csr.csr_matrix):
            XIST_logn = XIST_logn.toarray()
        data.obs['XIST_logn'] = XIST_logn[:,0]

        ygenes = data.var.index.isin(ygenes)
        yexpr = data.X[:,ygenes]
        if isinstance(yexpr, scipy.sparse.csr.csr_matrix):
            yexpr = yexpr.toarray()

    yexpr = yexpr.mean(axis = 1)
    data.obs['Ygene_logn'] = yexpr

    #plt.scatter(data.obs.XIST_logn, data.obs.Ygene_logn, alpha = 0.05)
    #plt.show()

    data.obs['XIST_bin'] = data.obs['XIST_logn'] > 0
    data.obs['Ygene_bin'] = data.obs['Ygene_logn'] > 0
    temp = pd.crosstab(index = data.obs.XIST_bin, columns = data.obs.Ygene_bin)
    #sns.heatmap(temp, annot = True)
    print(temp)
    
    data.obs['sex'] = 'UKN'
    data.obs.loc[data.obs.XIST_bin & ~data.obs.Ygene_bin , data.obs.columns == 'sex'] = 'F'
    data.obs.loc[data.obs.XIST_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'MF'
    data.obs.loc[~data.obs.XIST_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'M'
    
