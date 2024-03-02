from typing import Optional, Sequence, Union

import pandas as pd
from matplotlib.axes import Axes

from trackc.pa import trackcl_11
from trackc.tl._getRegionsCmat import GenomeRegion


def zoomin(
    ax: Optional[Axes] = None,
    raw_regions: Union[Sequence[str], str, None] = None,
    zoomin_regions: Union[Sequence[str], str, None] = None,
    colors: Union[Sequence[str], None] = None,
    alpha: float = 1,
    line_on: bool = False,
    fill: bool = True,
):
    """
    Plot zoomin track, support for multiple or reverse genome regions.

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    raw_regions: `str` | `str list`
        The raw genome regions, some of these regions will be selected to zoom in.
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``

        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), which means the reverse region
    zoomin_regions: `str` | `str list`
        regions to be zoomin, The format is the same as `raw_regions`.
        the regions of `zoomin_regions` should be located in the `raw_regions`
    colors: `str list`
        the colors of the zoomin_regions
    alpha: `float`
    line_on: `bool`
        Whether to display the line of zoomin plots
    fill: `bool`
        Whether to fill the zoomin plots

    Example
    -------
    >>> full_regions = "chr18:45000000-78077248"
    >>> zoom_regions = ['chr18:47400000-48280000', 'chr18:75280000-74030000']
    >>> neo_domain_regions = ['chr18:47950000-48280000', 'chr18:75280000-74850000']

    >>> ten = tc.tenon(figsize=(8,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=1, hspace=0.1)

    >>> tc.pl.zoomin(ax=ten.axs(0), raw_regions=full_regions, zoomin_regions=zoom_regions, line_on=True, fill=False, alpha=0.5)
    >>> tc.pl.zoomin(ax=ten.axs(1), raw_regions=zoom_regions, zoomin_regions=neo_domain_regions, line_on=False, fill=True, alpha=0.5)
    >>> tc.savefig('trackc_zoomin_track.pdf')
    """

    if isinstance(raw_regions, str):
        raw_regions = [raw_regions]
    if isinstance(zoomin_regions, str):
        zoomin_regions = [zoomin_regions]

    row_GRs = _region_pos(raw_regions)
    zoomin_GRs = _region_pos(zoomin_regions)
    zoomin_GRs["ors"] = None
    zoomin_GRs["ore"] = None
    zoomin_GRs_cp = zoomin_GRs.copy()
    for i, row in zoomin_GRs_cp.iterrows():
        start, end = _get_zoom_origin_pos(
            row["chrom"], row["fetch_start"], row["fetch_end"], row_GRs
        )
        ss = start
        ee = end
        if row["isReverse"] == True:
            ss = end
            ee = start
        zoomin_GRs.loc[i, "ors"] = ss
        zoomin_GRs.loc[i, "ore"] = ee

    zoomin_GRs = zoomin_GRs.reset_index()

    if colors == None:
        colors = trackcl_11


    if len(colors) < zoomin_GRs.shape[0]:
    	colors = colors * (int(zoomin_GRs.shape[0]/len(colors)) + 1)
    	print(len(colors))
    for i, row in zoomin_GRs.iterrows():
        x = row[["ps", "pe", "ore", "ors"]]
        y = [0, 0, 1, 1]

        if fill == True:
            ax.fill(x, y, color=colors[i], alpha=alpha)

        if line_on == True:
            ax.plot(
                (x[0], x[3]),
                (y[0], y[3]),
                color=colors[i],
                alpha=alpha,
                solid_capstyle="butt",
            )
            ax.plot(
                (x[1], x[2]),
                (y[1], y[2]),
                color=colors[i],
                alpha=alpha,
                solid_capstyle="butt",
            )

    ax.set_axis_off()
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])


def _region_pos(regions):
    GRs = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    GRs["len"] = GRs["fetch_end"] - GRs["fetch_start"]
    GRs["pos_e"] = GRs["len"].cumsum()
    GRs["pos_s"] = GRs["pos_e"] - GRs["len"]
    GRs_fulllen = GRs["len"].sum()
    GRs["ps"] = GRs["pos_s"] / GRs_fulllen
    GRs["pe"] = GRs["pos_e"] / GRs_fulllen
    return GRs


def _get_zoom_origin_pos(zoom_chrom, zooms, zoome, origin_pos_df):
    fulllen = origin_pos_df["len"].sum()
    origin_pos = origin_pos_df[origin_pos_df["chrom"] == zoom_chrom]
    origin_pos = origin_pos[
        (zooms >= origin_pos["fetch_start"]) & (zoome <= origin_pos["fetch_end"])
    ]

    start = None
    end = None
    if origin_pos.shape[0] == 0:
        print("{0}:{1}-{2} not in origion regions".format(zoom_chrom, zooms, zoome))
    else:
        ix = origin_pos.index[0]
        start_dis = origin_pos_df.loc[ix, "pos_s"] + abs(
            zooms - origin_pos_df.loc[ix, "start"]
        )
        end_dis = origin_pos_df.loc[ix, "pos_s"] + abs(
            zoome - origin_pos_df.loc[ix, "start"]
        )
        start = start_dis / fulllen
        end = end_dis / fulllen
    return start, end
