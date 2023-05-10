import matplotlib.pyplot as plt
from typing import List

def make_spec(width: float = 10,
              height: float = 3,
              height_ratios: List[float]=[1],
              hspace: float = 0,

              ):
    """\
    """
    
    nrows = len(height_ratios)
    if nrows <= 0:
        fig, ax = plt.subplots(1, 1, figsize=(width, height) )
        return fig, ax
    
    fig = plt.figure(figsize=(width, height))
    spec = fig.add_gridspec(nrows=nrows, ncols=1,
                           left=None, bottom=None, right=None, top=None,
                           wspace=0, hspace=hspace, 
                           width_ratios=None, 
                           height_ratios=height_ratios)
    
    axs = [None] * nrows
    for i, v in enumerate(height_ratios):
        axs[i] = fig.add_subplot(spec[i])

    return fig, axs

# fig, axs = make_spec(height_ratios=[1,1.5,2.4], hspace=0.1)

class lego:
    def __init__(self, width=6, height=1) -> None:
        pass
        fig, ax = plt.subplots(1, 1, figsize=(width, height) )
        self.fig = fig
        self.ax = ax
    def add(self, pos='top', height_ratio=1)



def savefig(outfile):
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')

