# Global import of functions
# from .globimport import *

from .describe import * # Describe an anndata
from .get_from_raw import * # copies .raw.X to X
from .is_outlier import * # Outlier



from . import pl as plot
from . import io as ios # add notebook to anndata
from . import tl as tl
from . import pp as pp
import scanpy_minion

import scanpy as sc
sc.external.minion =  scanpy_minion
