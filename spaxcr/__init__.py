from . import tools as tl
from . import plot as pl
from . import _utils
from . import process as pp
from . import extra as ex


import sys

sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["pp","tl","pl",'ex']})