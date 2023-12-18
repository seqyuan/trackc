from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.collections import PatchCollection
from matplotlib.colors import Colormap
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib.colors import LogNorm #,CenteredNorm, SymLogNorm,PowerNorm, Normalize  
import numpy as np
from trackc.tl._getRegionsCmat import GenomeRegion
from .links import _plot_loop_arc
from .bigwig import _make_multi_region_ax
from trackc.pa import fruitpunch

def bdgmat_track(ax: Optional[Axes] = None,
                 bed: Union[pd.DataFrame, str, None]  = None,
            regions: Union[Sequence[str], str, None] = None,
            style: Union[str, None] = 'heatmap',
            logdata: bool = False, 
            minrange: Optional[float] = None,
            maxrange: Optional[float] = None,
            color: Union[Sequence[str], None] = 'tab:silver',
            cmap: Union[Colormap, str, None] = None,
            #show_names: Union[bool, None] = False,
            alpha: Union[float, None] = 1,
            label: Union[str, None] = None,
            label_fontsize: Union[int, None] = 12,
            tick_fontsize: Optional[int] = 8,
            tick_fl: Optional[str] ='%0.2f',
            colorbar_width: float = 0.3,
            ):
    """\
    Plot bedGraph matrix track, support for multiple or reverse genome regions.
    bedgraph matrix file like bedgraph, columns from 4th are values
    bedgraph matrix should be sorted by chrom, start

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    bed: `pd.DataFrame` | `str`
        If ``bed`` if a filepath, the file should have headers, The first three columns are [chrom start end]
        bedGraph matrix format: [chrom start end scorename1 scorename2 ...]
    regions: `str` | `str list`
        The genome regions to plot
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means you want to get the reverse region
    style: `str`
        heatmap or line, default heatmap
    logdata: `bool` | `bool list`
        whether log the data before plotting
    minrange: `float` | `float list`
        the minimum range of values used to show
    maxrange: `float` | `float list` 
        the maximum range of values used to show
    color: `str` or `list`
        the color of the plot for style: line, if color is a list, the line color will be set by matrix columns order
    cmap: `str` | `matplotlib.colors.Colormap`
        the colormap of the plot for style: heatmap
    show_h
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    #axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    if isinstance(bed, str)==True:
        bed=pd.read_table(bed, sep="\t", header=0)

    if isinstance(color, list)==False:
        color = [color]
    score_columns = bed.shape[1]-3
    if len(color) < score_columns:
        repeat_times = (line_GenomeRegions.shape[0] + len(color) - 1) // len(color)
        color = (color * repeat_times)[:score_columns]

    bed['chrom'] = bed['chrom'].astype(str)

    df2plot = pd.DataFrame()
    bed = bed.fillna(0)
    for ix, row in line_GenomeRegions.iterrows(): 
        bed2plot = bed[(bed['chrom']==row['chrom']) & (bed['end']>=row['fetch_start']) & (bed['start']<=row['fetch_end'])]
        if row['isReverse']==True:
            bed2plot = bed2plot.iloc[::-1]
        if df2plot.shape[0]==0:
            df2plot = bed2plot
        else:
            df2plot = pd.concat([df2plot, bed2plot], axis=0)

    norm = None
    if logdata:
        norm=LogNorm()
    clim = (minrange, maxrange)

    df2plot = df2plot.iloc[:,3:].T
    if style=="heatmap":
        caxes = ax.pcolormesh(df2plot, cmap=cmap, edgecolor='none', norm=norm, snap=True, linewidth=.001)
        _colorbar(ax, caxes)
        ax.set_xticks([])  # 隐藏刻度
        ax.set_xticklabels([])  # 隐藏刻度标签
        #if
        ax.set_yticks([i+0.5 for i in range(0, df2plot.shape[0])])  # 隐藏刻度
        ax.set_yticklabels(df2plot.index)  # 隐藏刻度标签
        
def _colorbar(axm,im, width="3%", height="50%"):
    axins = inset_axes(axm,
                   width=width,  # width = 5% of parent_bbox width
                   height=height,  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.01, 0.5, 1, 1),
                   bbox_transform=axm.transAxes,
                   borderpad=0,
                   )
    cbar=plt.colorbar(im, cax=axins, orientation='vertical')
    cbar.set_ticks([im.get_array().min(), im.get_array().max()])
    cbar.ax.yaxis.tick_right()
    cbar.ax.spines['top'].set_color('none')
    cbar.ax.spines['right'].set_color('none')
    cbar.ax.spines['bottom'].set_color('none')
    cbar.ax.spines['left'].set_color('none')

