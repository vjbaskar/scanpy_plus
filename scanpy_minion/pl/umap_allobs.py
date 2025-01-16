from ..globimport import *

def umap_allobs(data, wanted_obs, cmap, unwanted_obs = ['predicted_doublet', 'xist_bin', "Ygene_bin"], ncols = 4, plotw = 7, ploth = 6):
    """
    Plot umap of all obs
    If wanted_obs is not defined: Plot all obs sans unwanted_obs
    
    wanted_obs: list. list of obs column names
    cmap: matplotlib.colors.LinearSegmentedColormap
    unwanted_obs: Some obs like True,False will fail. Use `unwanted_obs` to remove them out
    ncols: number of cols
    plotw: plot size for a single subplot
    ploth: plot height for a single subplot
    
    Runs sc.pl.umap()
    returns none
    """

    #vals to plot
    obsvars = data.obs.columns

    if wanted_obs:
        obsvars = wanted_obs
    else:
        # Remove unwanted vars
        unwanted = ['predicted_doublet', 'xist_bin', "Ygene_bin"]
        obsvars = [ i  for i in obsvars if i not in unwanted ]
    
    
    # Setting ncols and calc nrows accordingly
    ncols = 4
    nrows = round(len(obsvars)/ncols) + 1

    # setting total figsizes
    plotw = plotw * ncols
    ploth = ploth * nrows


    #fig, ax = plt.subplots(nrows,ncols, figsize=(plotw,ploth))

    plt.figure(figsize=(plotw, ploth))
    for n, i in enumerate(obsvars):
        print(i, end = ' ')
        ax = plt.subplot( nrows, ncols, n + 1)
        sc.pl.umap(data, color = i, title = i, ax = ax, show = False, legend_loc="", cmap = cmap)
