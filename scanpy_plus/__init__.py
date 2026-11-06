# Global import of functions
# from .globimport import *

# scanpy-plus has no hard dependencies (light install). Fail early with a
# helpful message if the scanpy stack is missing.
try:
    import scanpy as sc
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "scanpy-plus requires scanpy. Install it yourself, or use: "
        "pip install 'scanpy-plus[scanpy]'  (or '[all]' for every optional tool)"
    ) from e

from .describe import * # Describe an anndata
from .get_from_raw import * # copies .raw.X to X
from .is_outlier import * # Outlier


from . import pp as pp
from . import pl as plot
from . import io as ios # add notebook to anndata
from . import tl as tl
from . import recipe as recipe


import scanpy_plus

sc.external.scanpy_plus = scanpy_plus
