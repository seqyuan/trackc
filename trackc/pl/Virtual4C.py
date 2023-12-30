from typing import Optional, Sequence, Union

import cooler
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.colors import Colormap

from trackc.pl.mapc import getData2Map, mapC
from trackc.tl._getRegionsCmat import extractContactRegions

fruitpunch = sns.blend_palette(["white", "red"], as_cmap=True)


def _plot4C_line_bar(
    data,
    ax,
    color,
    track_type,
    logdata,
    trim_range,
    minrange,
    maxrange,
):
    data, maxrange, minrange = getData2Map(
        data, maxrange=maxrange, minrange=minrange, trim_range=trim_range, inplace=False
    )
    print("maxrange:", maxrange, "minrange:", minrange)

    if logdata:
        data = np.log2(data)

    if track_type == "line":
        ax.plot(range(len(data[0, :])), data[0, :], color=color, solid_capstyle="butt")

    if track_type == "bar":
        ax.bar(
            x=range(len(data[0, :])),
            height=data[0, :],
            bottom=minrange,
            width=1,
            color=color,
            align="edge",
        )
    return (data, maxrange, minrange)


def _get_pets(df, binsize=10000):
    pets = df.reset_index()
    pets_df = pd.DataFrame(index=range(pets["cbins"].sum()))
    pets_df["chrom"] = "None"
    pets_df["bin_N"] = 0

    chroms = []
    bins = []

    for i, row in pets.iterrows():
        chroms.extend([row["chrom"]] * row["cbins"])

        bin_s = int(row["fetch_start"] / binsize)
        # if row['fetch_start']%binsize > 0:
        #    bin_s = bin_s

        bin_e = int(row["fetch_end"] / binsize)
        if row["fetch_end"] % binsize > 0:
            bin_e = bin_e + 1

        chrbins = list(range(bin_s, bin_e))

        if row["isReverse"] == True:
            chrbins = chrbins[::-1]
        bins.extend(chrbins)

    pets_df["chrom"] = chroms
    pets_df["chrbin"] = bins
    return pets_df


def virtual4C(
    clr: Union[cooler.Cooler, str],
    ax: Optional[Axes] = None,
    balance: bool = False,
    # divisive_weights = None,
    target: Union[str, None] = None,
    regions: Union[Sequence[str], str, None] = None,
    track_type: Union[str, None] = "line",
    color: Union[str, None] = "tab:blue",
    cmap: Union[Sequence[Colormap], str, None] = fruitpunch,
    target_color: Union[str, None] = None,
    target_name: Union[str, None] = None,
    logdata: bool = False,
    trim_range: float = 0.98,
    minrange: float = None,
    maxrange: float = None,
    label: Optional[str] = None,
    label_rotation: Union[int, None] = 0,
    label_fontsize: Optional[int] = 12,
):
    """
    Plot virtual4C track, support for multiple or reverse genome regions.

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    clr: `cooler.Cooler` | `str`
        cool format Hi-C matrix (https://github.com/open2c/cooler) or cool file path
        `cooler.Cooler` e.g. GM12878=cooler.Cooler('./GM12878.chr18.mcool::/resolutions/50000')
        `str` e.g. 'GM12878.chr18.mcool::/resolutions/50000' or 'GM12878.chr18.cool'
    balance: `bool`
        The ``'balance'`` parameters of ``coolMat.matrix(balance=False).fetch('chr6:119940450-123940450')``
    target: `str`
        Coordinates of the viewpoint. e.g. 'chr8:127735434-127735435'
    regions: `str` | `str list`
        The genome regions, which contact to the viewpoint
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``),
        which means you want to get the reverse region contacts
    track_type: `str`
        virtual4C types, you can choose one of ['line', 'bar', 'heatmap']
    color: `str`
        `line` or `bar` type virtual4C color
    cmap: `str` | `matplotlib.colors.Colormap`
        if track_type set `heatmap`, the cmap
    target_color: `str`
        virtual4C viewpoint color
    target_name: `str`
        viewpoint label
    logdata: `bool`
        do you want to log the data before plotting
    trim_range: `float`
        remove the extreme values by trimming the counts.[0,1]
    minrange: `float`
        the minimum range of values used to define the ylim or color palette
    maxrange: `float`
        the maximum range of values used to define the ylim or color palette
    label: `str`
        the title of the track, will show on the left
    label_rotation: `int`
        the label text rotation
    label_fontsize: `int`
        the label text fontsize

    Example
    -------
    >>> import trackc as tc
    >>> regions = ['chr8:127000000-129200000', 'chr14:96500000-99300000']
    >>> MYC_TSS = 'chr8:127735434-127735435'
    >>> ten = tc.tenon(figsize=(8,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> AML_1360 = cooler.Cooler('./GSM4604287_1360.iced.mcool::/resolutions/10000')
    >>> tc.pl.virtual4C(ax=ten.axs(0), clr=AML_1360, target=MYC_TSS, regions=regions,
                track_type='line', label='Virtual 4C', target_color='r')
    >>> tc.savefig('trackc_virtual4c.pdf')
    """
    if isinstance(clr, str):
        clr = cooler.Cooler(clr)

    data = extractContactRegions(
        clr, balance=balance, row_regions=target, col_regions=regions
    )
    cols = _get_pets(data.col_regions, clr.binsize)
    rows = _get_pets(data.row_regions, clr.binsize)
    cols["target"] = 0
    cols.loc[cols["chrbin"].isin(rows["chrbin"]), "target"] = 1
    targets_point = cols.query("target==1")

    target_x = 0
    target_height = 0
    target_bottom = 0

    if track_type in ["line", "bar"]:
        _, maxrange, minrange = _plot4C_line_bar(
            data=data.cmat,
            ax=ax,
            color=color,
            track_type=track_type,
            logdata=logdata,
            trim_range=trim_range,
            minrange=minrange,
            maxrange=maxrange,
        )

        # ax.bar_label(target_bar, label_type='edge')
        ax.text(
            targets_point.index[0],
            minrange + (maxrange - minrange) / 2,
            target_name,
            va="center",
            ha="right",
            rotation=90,
        )
        ax.set_ylim(minrange, maxrange)
        ax.set_xlim(0, data.cmat.shape[1] - 1)

        target_x = targets_point.index[0]
        target_height = maxrange
        target_bottom = minrange

    if track_type == "heatmap":
        mapC(
            ax=ax,
            mat=data.cmat,
            cmap=cmap,
            logdata=logdata,
            symmetric=True,
            trim_range=trim_range,
            minrange=minrange,
            maxrange=maxrange,
            map_type="square",
            ax_on=True,
            aspect="auto",
        )
        ax.set_aspect(aspect="auto")

        target_bottom = -0.5
        target_height = 1
        target_x = targets_point.index[0] - 0.5

    target_bar = ax.bar(
        x=target_x,
        height=target_height,
        bottom=target_bottom,
        width=1,
        color=target_color,
        align="edge",
        label=target_name,
    )

    ax.tick_params(bottom=False, top=False, left=False, right=True)
    ax.yaxis.tick_right()
    # ax.set_xticklabels("")
    ax.set_xticklabels("")
    ax.set_ylabel(
        label, fontsize=label_fontsize, rotation=label_rotation, ha="right", va="center"
    )
