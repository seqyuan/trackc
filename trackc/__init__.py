# the actual API
import warnings

warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

from trackc.gs import make_spec, savefig, tenon

__all__ = ["tl", "pl", "pa", "make_spec", "tenon"]
