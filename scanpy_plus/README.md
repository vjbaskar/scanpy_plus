# Scanpy misc functions

# Installation

```
pip install -U git+https://github.com/vjbaskar/scvm.git
```
Depends only on setuptools and plotly. It uses packages in the `scverse` but I have not hard `required` it in the package as it may overwrite your version of important packages.
## Usage
```
import scvm
```
## Help
```
help(scvm)
```

## Documentation
`adata.describe()`: Describes the type of data in adata.X = raw/norm/log/scaled

scvm.pl.umap3d(): 3D umap
