from ..globimport import *
import time

def umap_grid(data, color, nrows, ncols, figsize, **kwargs):
    """
    Plot a set of obs or genes as a grid
    returns: matplotlib fig obj
    kwargs passed to sc.pl.umap
    """
    fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize = figsize)
    axes = ax.ravel()
    for a in axes:
        a.set_axis_off()
    
    for i,c in enumerate(color):
        axes[i].set_axis_on()
        sc.pl.umap(data, color = c, ax = ax[i], show = False, **kwargs)
    fig.tight_layout()
