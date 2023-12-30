from typing import Sequence, Union

import bioframe as bf
from matplotlib.axes import Axes

from .zoomin import _region_pos


def vhighlight(
    axs: Union[Sequence[Axes], Axes, None] = None,
    regions: Union[Sequence[str], str, None] = None,
    light_regions: Union[Sequence[str], str, None] = None,
    colors: Union[Sequence[str], None] = "yellow",
    alpha: float = 0.3,
):
    """
    Plot vhighlight track, support for multiple or reverse genome regions.

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    regions: `str` | `str list`
        The raw genome regions, some of these regions will be selected to zoom in.
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``),
        which means the reverse region
    light_regions: `str` | `str list`
        regions to be zoomin, The format is the same as `regions`.
        the regions of `light_regions` should be located in the `regions`
    colors: `str list`
        the colors of the light_regions
    alpha: `float`

    Example
    -------
    >>> regions = ['18:47400000-48280000', '18:75280000-74030000']
    >>> light_regions = ['18:47950000-48280000', '18:75280000-74850000']

    >>> ten = tc.tenon(figsize=(6,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=1, hspace=0.1)

    >>> tc.pl.vhighlight(axs=[ten.axs(0), ten.axs(1)], colors=['y', 'b'], regions=regions, light_regions=light_regions)
    >>> tc.pl.multi_scale_track(ten.axs(2), regions=regions, scale_adjust='Mb', intervals=2)
    """

    if isinstance(regions, str):
        regions = [regions]
    if isinstance(light_regions, str):
        light_regions = [light_regions]

    row_GRs = _region_pos(regions)
    light_GRs = _region_pos(light_regions)
    len_sum = sum(row_GRs["len"])

    row_GRs = row_GRs[
        ["chrom", "fetch_start", "fetch_end", "isReverse", "len", "ps", "pe"]
    ]
    row_GRs.columns = ["chrom", "start", "end", "isReverse", "len", "ps", "pe"]

    light_GRs["light_width"] = light_GRs["len"] / len_sum
    light_GRs = light_GRs[["chrom", "fetch_start", "fetch_end", "light_width"]]
    light_GRs.columns = ["chrom", "start", "end", "light_width"]

    row_GRs = row_GRs.reset_index()
    light_GRs = light_GRs.reset_index()
    del light_GRs["index"]
    del row_GRs["index"]

    light_r = bf.overlap(light_GRs, row_GRs)
    light_r = light_r.query("chrom_==chrom_")

    light_r["light_s"] = 0
    light_r_cp = light_r.copy()

    for i, row in light_r_cp.iterrows():
        if row["isReverse_"] == True:
            light_r.loc[i, "light_s"] = (
                row["ps_"] + abs(row["end"] - row["end_"]) / len_sum
            )
        else:
            light_r.loc[i, "light_s"] = (
                row["ps_"] + abs(row["start"] - row["start_"]) / len_sum
            )

    if isinstance(colors, list) == False:
        colors = [colors]

    if len(colors) < light_r.shape[0]:
        repeat_times = (light_r.shape[0] + len(colors) - 1) // len(colors)
        colors = (colors * repeat_times)[: light_r.shape[0]]

    for ix, aix in enumerate(axs):
        # trans = transforms.blended_transform_factory(
        #    aix.transAxes, aix.transAxes)
        aix.bar(
            x=light_r["light_s"],
            height=1,
            width=light_r["light_width"],
            transform=aix.transAxes,
            align="edge",
            color=colors,
            # zorder=50,
            clip_on=False,
            alpha=alpha,
        )
