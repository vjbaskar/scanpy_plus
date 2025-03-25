from ..globimport import *

"""
Allows to create a persistent png image from
the last plotted matplotlib plot widget
Before capturing do %matplotlib widget
"""
def perma_plot():
    import base64
    from io import BytesIO
    from IPython.display import HTML
    # TODO: maybe pass parameters for savefig
    #       to control quality/type of img
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    encoded_string = base64.b64encode(buffer.read()).decode('utf-8')
    html_string = '<img src=\'data:image/png;base64,{}\'>'.format(encoded_string)
    return HTML(html_string)
