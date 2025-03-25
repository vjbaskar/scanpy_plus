from ..globimport import *

def plotly_umap(data, color_col, obs_cols = None, hover_data_cols = None,  write_to_file = 'umap.html', point_size = 3, **kwargs):
    
    """
    Generates a plotly umap
    
    Arguments
    ---------
    data: adata or mdata
    color_col: str, data.obs column name 
    obs_cols: [list], data.obs column names to be added to the umap plot
    hover_data_cols: [list], subset of obs_cols that should appear in mouse hover mode
    write_to_file: str, write to html file
    kwargs: Passed to the plotly.express.scatter function
    
    Example:
    # Get colors for clusters from the mdata
    my_colors = dict(zip(mdata.obs.annotation.cat.categories, mdata.uns['annotation_colors']))
    my_colors # a dict {'3-CD4_T': '#ff4a46'}
    
    fig = plotly_umap(mdata, color_col = 'annotation', obs_cols = ['annotation', 'celltype'], hover_data_cols = ['annotation', 'celltype'], width=1100, height=900, template="simple_white", color_discrete_map = my_colors)
    
    """
    
    import plotly.express as px
    from IPython.display import HTML

    df = pd.DataFrame(data.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=data.obs.index)
    try:
        if type(color_col) is str:
            obs_cols.append(color_col)
        elif type(color_col) is list:
            obs_cols.extend(color_col)
    except NameError:
        obs_cols = color_col
    
    try:
        obs_cols.extend(hover_data_cols)
    except NameError:
        hover_data_cols = color_col
    
    obs_cols = list(set(obs_cols))  
    df = df.join(data.obs[obs_cols])
    fig = px.scatter(df, x="UMAP1", y="UMAP2", color=color_col,
                     hover_data=hover_data_cols,  **kwargs)
    fig.update_traces(marker=dict(size=point_size), selector=dict(mode='markers'))
    fig.update_layout(legend= {'itemsizing': 'constant'}) # Keep figure legend marker size constant not according to point size
    fig.show()
    fig.write_html(write_to_file)
    return fig