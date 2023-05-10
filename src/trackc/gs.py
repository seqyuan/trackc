import matplotlib.pyplot as plt
from typing import List
from matplotlib.axes import Axes

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

class _lego_ax:
    bottom = 0
    top = 0
    track_height = 0
    hspace = 0
    ax = None
    def __init__(self, ax_ref: Axes, bottom: float, track_height: float, hspace: float):
        self.ax = ax_ref.inset_axes([0, bottom, 1, track_height])
        self.bottom = bottom
        self.top = bottom + hspace
        self.track_height = track_height
        self.hspace = hspace
    

class lego:
    lego_axs = []
    fig = None
    ax = None
    def __init__(self, width=6, height=1):
        fig, ax = plt.subplots(1, 1, figsize=(width, height) )
        ax.set_axis_off()
        self.fig = fig
        self.ax = ax
    def show(self):
       self.fig.show()

    def add(self, pos='top', height_ratio=1, hspace=0):
        ref_lego_ax = _lego_ax(self.ax, 0, 0, 0)
        if len(self.lego_axs) != 0:
            ix = -1
            if pos=="top":
                ix = 0 
            ref_lego_ax = self.lego_axs[ix]
        else:
            hspace = 0
        
        bottom = ref_lego_ax.bottom - hspace - height_ratio
        if pos == 'top':
            bottom = ref_lego_ax.top + hspace

        lego_ax_demo = _lego_ax(self.ax, bottom, height_ratio, hspace)
        if pos == 'top':
            self.lego_axs.insert(0, lego_ax_demo)
        else:
            self.lego_axs.append(lego_ax_demo)

    def remove(self, pos='top'):
        if len(self.lego_axs) >=1:
            if pos == 'top':
                self.lego_axs.remove(self.lego_axs[0])
            elif pos == 'bottom':
                self.lego_axs.remove(self.lego_axs[-1])
            else:
                pass


def savefig(outfile):
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')

