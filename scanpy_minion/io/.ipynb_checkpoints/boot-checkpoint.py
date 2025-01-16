import anndata

def push(adata, notebook_name, save_filename):
    """
    Push a notebook into the adata
    notebook_name: str. 
    """
    from datetime import datetime
    import os
    import subprocess

    notebook_basename=os.path.splitext(os.path.basename(notebook_name))[0]
    #cmd = ['/home/jovyan/my-conda-envs/singlecell/bin/jupytext', '--to', 'py:percent', notebook_name]
    cmd = ['/usr/local/bin/singularity', 'exec', 
           
           '--bind', os.getcwd()+':/mnt',
           '/nfs/team298/vm11/soft/__singularity/jupytext.sif',
           'jupytext',
           '--to', 'py:percent',
           '/mnt/'+notebook_name
          ]
    subprocess.run(cmd, check=True)
    notebook_py = notebook_basename + ".py"
    with open(notebook_py, 'r') as f:
        x = f.readlines()
    adata.uns['code_base'] = dict()
    adata.uns['code_base']['code'] = x
    adata.uns['code_base']['notebook_basename'] = notebook_basename
    adata.uns['code_base']['time_of_writing'] = datetime.now().strftime("%Y-%m-%d:%H:%M:%S")
    adata.write(save_filename)

    
def pull(adata, save_notebook = None):
    """
    Pull notebook from data
    You should use push first
    Args:
    adata: anndata self
    save_notebook: str. file_suffix to save file
                    default = None
    """
    
    import os
    x = adata.uns['code_base']['code']
    if save_notebook == None:
        notebook_basename = adata.uns['code_base']['notebook_basename']
    else:
        notebook_basename = save_notebook
    py_name = notebook_basename + ".py"
    notebook_name = notebook_basename + "_reconstructed.ipynb"
    print(f"Writing to {notebook_name}")
    with open(py_name, "w") as f:
        for i in x:
            print(i, file=f)
    cmd = ['/usr/local/bin/singularity', 'exec', 
       '--bind', os.getcwd()+':/mnt',
       '/nfs/team298/vm11/soft/__singularity/jupytext.sif',
       'jupytext',
           '--from', 'py:percent',
       '--to', 'ipynb',
           '-o', '/mnt/' + notebook_name, 
       '/mnt/'+py_name
      ]
        
    os.system(" ".join(cmd))
        #subprocess.run(cmd,  check=True)
anndata.AnnData.push = push
anndata.AnnData.pull = pull