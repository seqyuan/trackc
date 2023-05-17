import matplotlib.pyplot as plt
from typing import List, Union
from matplotlib.axes import Axes

def make_spec(width: float = 10,
              height: float = 3,
              height_ratios: Union[List[float], None] = None,
              width_ratios: Union[List[float], None] = None,
              hspace: float = 0,
              wspace: float = 0,

              ):
    """\
    """
    
    nfig = 0
    nrows = 1
    ncols = 1
    ratios = height_ratios
    if height_ratios != None and isinstance(height_ratios, list):
        nfig = len(height_ratios)
        width_ratios = None
        nrows = nfig

    if width_ratios != None and isinstance(width_ratios, list):
        nfig = len(width_ratios)
        height_ratios = None
        ncols = nfig
        ratios = width_ratios

    if nfig <= 0:
        fig, ax = plt.subplots(1, 1, figsize=(width, height) )
        return fig, ax
    
    fig = plt.figure(figsize=(width, height))
    spec = fig.add_gridspec(nrows=nrows, ncols=ncols,
                           left=None, bottom=None, right=None, top=None,
                           wspace=wspace, hspace=hspace, 
                           width_ratios=width_ratios, 
                           height_ratios=height_ratios)
    
    axs = [None] * nfig
    for i, v in enumerate(ratios):
        axs[i] = fig.add_subplot(spec[i])

    return fig, axs

# fig, axs = make_spec(height_ratios=[1,1.5,2.4], hspace=0.1)

class mortise:
    bottom = 0
    top = 0
    track_height = 0
    hspace = 0
    ax = None
    def __init__(self, ax_ref: Axes, bottom: float, track_height: float, hspace: float):
        self.ax = ax_ref.inset_axes([0, bottom, 1, track_height])
        self.bottom = bottom
        self.top = bottom + track_height
        self.track_height = track_height
        self.hspace = hspace
    
class tenon:
    mortises = []
    fig = None
    ax = None
    def __init__(self, width=7, height=1 ):
        fig, ax = plt.subplots(1, 1, figsize=(width, height) )
        ax.set_axis_off()
        self.fig = fig
        self.ax = ax
        
    def show(self):
       for i, v in enumerate(self.mortises):
           #v.ax.set_ylabel('mortises[{0}]'.format(i), rotation=0)
           v.ax.text(0.5, 0.5, '.axs({0})'.format(i), ha='center', va='center')
           v.ax.set_xticklabels([])
           v.ax.set_xticks([])
       self.fig.show()

    def add(self, pos='top', height=1, hspace=0):
        if len(self.mortises) != 0:
            ix = -1
            if pos=="top":
                ix = 0 
            ref_lego_ax = self.mortises[ix]
        else:
            ref_lego_ax = mortise(self.ax, 0, height, 0)
            self.mortises = [ref_lego_ax]
            return
        
        bottom = ref_lego_ax.bottom - hspace - height
        if pos == 'top':
            bottom = ref_lego_ax.top + hspace

        lego_ax_demo = mortise(self.ax, bottom, height, hspace)
        if pos == 'top':
            self.mortises.insert(0, lego_ax_demo)
        else:
            self.mortises.append(lego_ax_demo)

    def remove(self, pos='top'):
        if len(self.mortises) >=1:
            if pos == 'top':
                self.mortises.remove(self.mortises[0])
            elif pos == 'bottom':
                self.mortises.remove(self.mortises[-1])
            else:
                pass
    def axs(self, index):
        return self.mortises[index].ax

def savefig(outfile):
    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight')

