# Global import of functions
# from .globimport import *

from .describe import * # Describe an anndata
from .get_from_raw import * # copies .raw.X to X
from .is_outlier import * # Outlier


from . import pp as pp
from . import pl as plot
from . import io as ios # add notebook to anndata
from . import tl as tl
from . import recipe as recipe


import scanpy_plus

import scanpy as sc
sc.external.scanpy_plus =  scanpy_plus
