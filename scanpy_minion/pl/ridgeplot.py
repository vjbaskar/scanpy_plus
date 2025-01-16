from ..globimport import *

def ridgeplot(adata, x = 'log1p_total_counts', split_by = 'pool', title = None, alpha = 1, height = 1.2, aspect = 9, overlap = 0.3):
    
    """
    Generates ridgeplot akin to Seurat.
    adata: anndata obj
    x: Which obs column with continuous variable to plot
    split_by: Which obs column to split by. For eg. Leiden clusters
    title: str
    alpha: transparency of the kdeplot
    height: height of the plot
    aspect: aspect
    overlap: overlap between two kdeplots.
    returns: None
    """
    
    if title == None:
        title = x + " | " + split_by
    
    obsdata = adata.obs
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
    g= sns.FacetGrid(obsdata, row = split_by, hue = split_by, aspect=aspect, height=height)
    g.map_dataframe(sns.kdeplot, x=x, fill=True, alpha=alpha)
    g.map_dataframe(sns.kdeplot, x=x, color='black')
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, color='black', fontsize=13,
                ha="left", va="center", transform=ax.transAxes)
    g.map(label, split_by)
    #g.fig.subplots_adjust(hspace=-.3)
    g.fig.subplots_adjust(hspace=-overlap)

    g.set_titles("")
    g.set_ylabels("")
    g.set(yticks=[])
    g.despine(left=True)
    plt.suptitle(title, y=0.98)

