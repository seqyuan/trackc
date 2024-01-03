from typing import Any, Callable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
import pyBigWig
from trackc._utils import LOGGER
from trackc.tl._getRegionsCmat import GenomeRegion
import warnings
warnings.filterwarnings("ignore")

def _make_multi_region_ax(ax, lineGenomeRegions):
    lineGenomeRegions["len"] = (
        lineGenomeRegions["fetch_end"] - lineGenomeRegions["fetch_start"]
    )
    lineGenomeRegions["ax_ratio"] = (
        lineGenomeRegions["len"] / lineGenomeRegions["len"].sum()
    )
    lineGenomeRegions["ax_x"] = (
        lineGenomeRegions["ax_ratio"].cumsum(axis=0) - lineGenomeRegions["ax_ratio"]
    )
    axs = [
        ax.inset_axes([row["ax_x"], 0, row["ax_ratio"], 1])
        for i, row in lineGenomeRegions.iterrows()
    ]
    for axi in axs:
        axi.axis("off")

    return axs


def bw_track(
    bw,
    ax: Optional[Axes] = None,
    regions: Union[Sequence[str], str, None] = None,
    binsize: Optional[int] = 50000,
    style: Optional[str] = "bar",
    summary_type: Union[str, None] = "mean",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    primary_col: Union[Sequence[str], None] = "#3271B2",
    secondary_col: Union[Sequence[str], None] = "#FBD23C",
    alpha: Optional[float] = 1.0,
    invert_y: Optional[bool] = False,
    label: Optional[str] = None,
    label_rotation: Union[int, None] = 0,
    label_fontsize: Optional[int] = 9,
    tick_fontsize: Optional[int] = 7,
    tick_fl: Optional[str] = "%0.2f",
    ax_on: bool = False,
):
    """
    Plot bigwig signal track, support for multiple or reverse genome regions.

    Parameters
    ----------
    bw: `pyBigWig.open` query object, or bigwig file path
    ax: :class:`matplotlib.axes.Axes` object
    regions: `str` | `str list`
        The genome regions to show the signal.
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``),
        which means the reverse region
    binsize: `int`
        binsize divided to computing signal summary statistics
    type: `str`
        plot type, default='bar', options=['bar', 'line']
    summary_type: `str`
        Summary type (mean, min, max, coverage, std), default 'mean'.
    vmin: `float`
        the minimum range of values used to define the ylim
    vmax: `float`
        the maximum range of values used to define the ylim
    primary_col: `str`
        the signal bar color
    secondary_col: `str`
        the signal bar color for negative values
    alpha: `float`
        alpha of plot color
    invert_y: `bool`
        whether reverse the y-axis
    label: `str`
        the title of the track, will show on the left
    label_rotation: `int`
        the label text rotation
    label_fontsize: `int`
        the label text fontsize
    tick_fontsize: `int`
        values range ticks text fontsize
    tick_fl: `str`
        values range ticks retains a few decimal places
    ax_on: `bool`
        whether show the spines

    Example
    -------
    >>> import trackc as tc
    >>> import pyBigWig
    >>> H3K27ac = pyBigWig.open('./GSM4604189.bigwig')

    >>> ten = tc.tenon(figsize=(8,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=1, hspace=0.2)

    >>> regions = ['chr8:127000000-129200000', 'chr14:96500000-99300000']
    >>> tc.pl.bw_track(H3K27ac, ten.axs(0), regions=regions, vmax=20, label='H3K27ac', binsize=10000, color='tab:blue')
    >>> tc.pl.bw_track(H3K27ac, ten.axs(1), regions=regions, vmax=5, label='H3K27ac', binsize=10000, invert_y=True, ax_on=True)
    >>> tc.savefig('trackc_bigwig_track.pdf')
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat(
            [GenomeRegion(i).GenomeRegion2df() for i in regions]
        )
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    if isinstance(bw, str) == True:
        bw = pyBigWig.open(bw)

    axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    if isinstance(primary_col, list) == False:
        primary_col = [primary_col]
    if len(primary_col) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(primary_col) - 1) // len(primary_col)
        primary_col = (primary_col * repeat_times)[: line_GenomeRegions.shape[0]]

    if isinstance(secondary_col, list) == False:
        secondary_col = [secondary_col]
    if len(secondary_col) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(secondary_col) - 1) // len(secondary_col)
        secondary_col = (secondary_col * repeat_times)[: line_GenomeRegions.shape[0]]

    min_y = 0
    max_y = 0
    plot_bottom_line = True

    for i, row in line_GenomeRegions.iterrows():
        bins = int(row["len"] / binsize)
        plot_list = bw.stats(
            row["chrom"],
            row["fetch_start"],
            row["fetch_end"],
            type=summary_type,
            nBins=bins,
        )
        plot_list = [0 if v is None else v for v in plot_list]
        if style == "line":
            axs[i].plot(range(0, bins), plot_list, color=primary_col[i], alpha=alpha)
        else:
            plotvalues = pd.DataFrame({'v':plot_list}, index=range(0, bins))
            pos_value = plotvalues.query('v>=0')
            neg_value = plotvalues.query('v<0')
            if neg_value.shape[0]>0:
                plot_bottom_line = False 

            axs[i].bar(
                x=pos_value.index,
                height=pos_value['v'],
                width=1,
                #bottom=[0] * (pos_value.shape[0]),
                bottom=0,
                color=primary_col[i],
                align="edge",
                edgecolor=None,
                alpha=alpha,
            )
            
            axs[i].bar(
                x=neg_value.index,
                height=neg_value['v'],
                width=1,
                #bottom=[0] * (neg_value.shape[0]),
                bottom=0,
                color=secondary_col[i],
                align="edge",
                edgecolor=None,
                alpha=alpha,
            )

        right, left = bins, 0
        if row["isReverse"] == True:
            left, right = bins, 0
        axs[i].set_xlim(left, right)

        if min_y > min(plot_list):
            min_y = min(plot_list)

        if max_y < max(plot_list):
            max_y = max(plot_list)

    if vmin == None:
        vmin = min_y
    if vmax == None:
        vmax = max_y

    for axi in axs:
        if invert_y:
            axi.set_ylim(vmax, vmin)
        else:
            axi.set_ylim(vmin, vmax)
    va = "top"
    if invert_y:
        va = "bottom"
        ax.set_ylim(vmax, vmin)
    else:
        ax.set_ylim(vmin, vmax)

    if plot_bottom_line:
        ax.text(
            0,
            vmax,
            " [{0}, {1}]".format(tick_fl % vmin, tick_fl % vmax),
            verticalalignment=va,
            fontsize=tick_fontsize,
        )
        ax.set_yticks([])
        ax.set_yticklabels("")

    else:
        ax.set_yticks([vmin, 0, vmax])
        ax.set_yticklabels([f'{tick_fl % vmin}', '0', f'{tick_fl % vmax}'], fontsize=tick_fontsize)
        if invert_y:
            ax.set_yticks([vmax, 0, vmin])
            ax.set_yticklabels([f'{tick_fl % vmax}', '0', f'{tick_fl % vmin}'], fontsize=tick_fontsize)
            
    ax.set_ylabel(
        label,
        fontsize=label_fontsize,
        rotation=label_rotation,
        horizontalalignment="right",
        verticalalignment="center",
    )

    if ax_on == False:
        spines = ["top", "bottom", "left", "right"]
        if invert_y == True:
            if plot_bottom_line==True:
                del spines[0]
            else:
                del spines[2]
        else:
            if plot_bottom_line==True:
                del spines[1]
            else:
                del spines[2]
        for i in spines:
            ax.spines[i].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels("")
    

def bw_compartment(
    compartment_bw,
    ax,
    chrom,
    start,
    end,
    ylabel,
    xticklabel=False,
    Acolor="#3271B2",
    Bcolor="#FBD23C",
    binsize=100000,
):
    chrsize = compartment_bw.chroms()[chrom]
    xbins = int(chrsize / binsize)
    if chrsize % binsize > 0:
        xbins = xbins + 1

    plot_list = compartment_bw.stats(chrom, 0, chrsize, nBins=xbins)
    mat = pd.DataFrame({"pc1": plot_list})
    mat["start"] = mat.index * binsize
    mat["width"] = binsize
    mat.loc[mat.index[-1], "length"] = chrsize % binsize
    mat.loc[mat[np.isnan(mat.pc1)].index, "pc1"] = 0

    plus = mat[mat["pc1"] > 0]
    minux = mat[mat["pc1"] <= 0]

    ax.bar(
        x=list(plus["start"]),
        height=plus["pc1"],
        width=plus["width"],
        bottom=[0] * (plus.shape[0]),
        color=Acolor,
        align="edge",
        edgecolor=Acolor,
        label="A",
    )
    ax.bar(
        x=list(minux["start"]),
        height=minux["pc1"],
        width=minux["width"],
        bottom=[0] * (minux.shape[0]),
        color=Bcolor,
        align="edge",
        edgecolor=Bcolor,
        label="B",
    )
    # ax.bar(0,height=0,color="#E27678",align="edge",edgecolor="#E27678",label="A2B")
    # ax.bar(0,height=0,color="#85AFBD",align="edge",edgecolor="#85AFBD",label="B2A")

    # ax.set_xlim(0, chrsize)
    ax.grid(False)
    ax.tick_params(bottom=False, top=False, left=True, right=False)  # 去掉tick线
    ax.spines["left"].set_color("k")
    ax.spines["left"].set_linewidth(1)
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.plot(
        [0, chrsize],
        [0, 0],
        "-",
        label="",
        linewidth=1,
        color="black",
        solid_capstyle="butt",
    )
    # ax.set_yticklabels('')
    ax.set_ylabel(
        ylabel,
        fontsize=10,
        rotation="horizontal",
        horizontalalignment="right",
        verticalalignment="center",
    )
    # ax.set_ylim([-1,1])
    # ax.set_yticks([-1, 1])
    ax.set_yticklabels([-1, 1], fontsize=8)
    if xticklabel == False:
        ax.set_xticklabels("")

    ax.set_xlim(start, end)
