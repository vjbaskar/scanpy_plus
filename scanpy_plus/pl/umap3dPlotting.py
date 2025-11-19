from ..globimport import *

def umap3dPlotting(adata, color = 'celltype', cols = ['leiden', 'condition'],
                   title = '3D UMAP', marker_size = 1, 
                   hovermode = True, output_html = 'outpt.html', figshow = False):
    """
    Plot umap 3d using plotly express.
    You have to run umap with 3d first.
    It expects the coordinates in anndata.obsm['X_umap'] which is the default.
    returns a data frame that it plots.
    
    Options:
    adata: anndata
    color: str. color umap by which col
    cols: list. cols to include in the dataframe
    title: str 
    marker_size: float. Size of points
    hovermode: boolean. Whether or not to have mouse over display
    output_html: str. 
    figshow: boolean. Show figure
    
    """
    import plotly.express as px
    print("Generating data frame")
    umap_matx = adata.obsm['X_umap'].copy()
    umap_df = pd.DataFrame(umap_matx)
    umap_df.index = adata.obs_names
    #umap_df[color] = adata.obs[color]
    if type(cols) == str:
        cols = [cols]

    for i in cols + [color]:
        print(i)
        umap_df[i] = adata.obs[i].copy()
    print(umap_df.head())
    umap_df.columns = ['X1','X2','X3'] + cols + [color]
    #temp = umap_df.head(10000)
    temp = umap_df
    print("Getting color info")
    temp2 = adata.uns[color + '_colors'] # Get colors from scanpy adata
    colordict = dict(zip( list(adata.obs[color].cat.categories), temp2)) # Create a dict of names: color
    
    print("Plotly express")
    
    #df = px.data.iris()
    fig = px.scatter_3d(temp, 
                        x='X1', y='X2', z='X3',
                        color=color, 
                        width=1000, height=800,
                       color_discrete_map=colordict
                       )
    #fig.update_traces(marker_line_color('rgb(0, 2, 1)'))
    fig.update_layout(title=title, autosize=False,
                      width=800, height=800)
    fig.update_traces(marker=dict(size=marker_size, 
                                  line=dict(width=0,
                                            color='gray')),
                      selector=dict(mode='markers'))
    fig.update_layout(template = 'plotly_white')
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.2
    ))
    #fig.update_layout(hovermode=hovermode)
    if figshow == True:
        print(figshow)
        fig.show()
    fig.write_html(output_html)
    return(temp)


