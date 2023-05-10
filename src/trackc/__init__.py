# the actual API
from . import tools as tl
from . import plotting as pl
from . import palettes as pa
from .gs import (
    make_spec,
    lego,
    savefig
)

import sys
sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['tl', 'pl', 'pa']})
from ._utils import annotate_doc_types
annotate_doc_types(sys.modules[__name__], 'trackc')
del sys, annotate_doc_types

