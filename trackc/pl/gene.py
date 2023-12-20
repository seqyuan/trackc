from typing import Optional, Sequence, Union

import pandas as pd
from matplotlib.axes import Axes

from trackc.tl._getRegionsCmat import GenomeRegion

from .bigwig import _make_multi_region_ax


def gene_track(
    ax: Optional[Axes] = None,
    bed12: Union[pd.DataFrame, str] = None,
    regions: Union[Sequence[str], str, None] = None,
    track_type: Union[str, None] = "gene",
    show_label: Union[bool, Sequence[str], str] = True,
    pos_strand_gene_color: Union[str, None] = "#3366CC",
    neg_strand_gene_color: Union[str, None] = "#EECFA1",
    line: Union[int, None] = 1,
    gene_fontszie: Union[int, None] = 7,
    label: Optional[str] = None,
    label_rotation: Union[int, None] = 0,
    label_fontsize: Optional[int] = 12,
    ax_on: bool = False,
):
    """
    Plot gene track, support for multiple or reverse genome regions.

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    bed12: `pd.DataFrame` | `str`
        gene annotation bed12 format `DataFrame` or `filepath`
        Bed12 files can be converted from GTF using `gtf2bed4trackc`.
        https://trackc.readthedocs.io/en/latest/tutorials/fileformats/bed12.html#gtf2bed4trackc
    regions: `str` | `str list`
        genome regions, format: `chrom:start-end`
        examples: ['chr18:47950000-48280000', 'chr18:75280000-74850000'] or "chr18:45000000-78077248"
        If the start is bigger than end, the genome region will be reversed
    track_type: `str`
        you can select one of the options: `gene` or `dendity`
        gene: gene track style
        dendity: gene density style. Under development
    show_label: `bool` | `str` | `str list`
        If the value is `False`, the gene name will not show
        If want show one gene, and hide others, just set the gene or gene list as the value, eg: `PIBF1` | `['PIBF1', 'KLF5']`
    pos_strand_gene_color: `str`
        positive strand gene name color
    neg_strand_gene_color: `str`
        negative strand gene name color
    line: `int`
        rows occupied by the genes in the region plotted
    gene_fontszie: `int`
        gene label fontszie
    label: `str`
        the title of the track, will show on the left
    label_rotation: `int`
        the label text rotation
    label_fontsize: `int`
        the label text fontsize
    ax_on: `bool`
        If True, top, left and right spines will show

    Returns
    -------
    None

    Example
    -------
    >>> import trackc as tc
    >>> regions = ['chr18:47950000-48280000', 'chr18:75280000-74850000']
    >>> gene_bed12 = '/path/GRCh38.84.bed12'

    >>> fig, axs = tc.make_spec(figsize=(7,2), height_ratios=[1])
    >>> tc.pl.gene_track(gene_bed12, ax=axs[0], regions=regions, line=12)
    >>> tc.savefig('trackc_gene_track.pdf')
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat(
            [GenomeRegion(i).GenomeRegion2df() for i in regions]
        )
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    ax.set_ylabel(
        label,
        fontsize=label_fontsize,
        rotation=label_rotation,
        horizontalalignment="right",
        verticalalignment="center",
    )

    ax.set_xticks([])
    ax.set_xticklabels("")
    ax.set_yticks([])
    ax.set_yticklabels("")
    if ax_on == False:
        spines = ["top", "bottom", "left", "right"]
        for i in spines:
            ax.spines[i].set_color("none")

    if isinstance(bed12, str):
        bed12 = pd.read_table(bed12, sep="\t", header=None)

    bed12 = bed12.iloc[:, 0:12]
    bed12.columns = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    bed12["blockSizes"] = bed12["blockSizes"].str.rstrip(",")
    bed12["blockStarts"] = bed12["blockStarts"].str.rstrip(",")
    bed12["chrom"] = bed12["chrom"].astype(str)

    for ix, row in line_GenomeRegions.iterrows():
        if track_type == "gene":
            _plot_gene(
                axs[ix],
                bed12,
                row["chrom"],
                row["fetch_start"],
                row["fetch_end"],
                needReverse=row["isReverse"],
                show_label=show_label,
                pos_strand_gene_color=pos_strand_gene_color,
                neg_strand_gene_color=neg_strand_gene_color,
                line=line,
                fontszie=gene_fontszie,
                ax_on=ax_on,
            )
        if track_type == "density":
            print("This gene type is developping")
        else:
            pass


def _plot_gene(
    ax,
    gene_bed,
    chrom,
    start,
    end,
    needReverse=False,
    show_label=True,
    pos_strand_gene_color="#3366CC",
    neg_strand_gene_color="#EECFA1",
    line=1,
    fontszie=5,
    ax_on=False,
):
    gene_bed = gene_bed[gene_bed["chrom"] == chrom]
    gene_bed_plot = gene_bed[
        ((gene_bed["start"] >= start) & (gene_bed["start"] <= end))
        | ((gene_bed["end"] >= start) & (gene_bed["end"] <= end))
    ]
    gene_bed_plot = gene_bed_plot.sort_values(by="end")
    # print(gene_bed_plot

    plot_gene_num = gene_bed_plot.shape[0]

    ii = 0
    head_length = (abs(end - start) / (line + 2)) / 5
    if line <= 3:
        head_length = (abs(end - start) / (line * 3)) / 10

    for i, row in gene_bed_plot.iterrows():
        # col = pos_strand_gene_color
        text_col = pos_strand_gene_color

        if row["strand"] == "-":
            # col = neg_strand_gene_color
            text_col = neg_strand_gene_color

        # text_col = col
        plot_y = ii % line

        ax.plot(
            (row["start"], row["end"]),
            (plot_y + 0.5, plot_y + 0.5),
            color="k",
            linewidth=1,
            solid_capstyle="butt",
        )
        starts = [int(x) for x in row["blockStarts"].split(",")]
        widths = [int(x) for x in row["blockSizes"].split(",")]

        ax.bar(
            x=starts,
            height=0.4,
            width=widths,
            bottom=plot_y + 0.3,
            edgecolor="k",
            linewidth=1,
            align="edge",
            color="k",
        )

        if row["start"] < start:
            row["start"] = start
        if row["end"] > end:
            row["end"] = end

        arrow_s = row["end"]
        dx = 0.3
        if row["strand"] == "-":
            arrow_s = row["start"]
            dx = -0.1
        ax.arrow(
            arrow_s,
            plot_y + 0.5,
            dx,
            0,
            overhang=0.5,
            width=0,
            head_width=0.28,
            head_length=head_length,
            length_includes_head=False,
            color=text_col,
            linewidth=0.5,
        )
        if isinstance(show_label, bool):
            if show_label == False:
                ii += 1
                continue
        if isinstance(show_label, str):
            if row["name"] != show_label:
                ii += 1
                continue
        if isinstance(show_label, list):
            if row["name"] not in show_label:
                ii += 1
                continue

        if (row["name"] in gene_bed_plot.iloc[-4:, :]["name"]) or (
            int(line / (ii + 1)) < 2
        ):
            ha = "right"
            genename = row["name"] + "  "
            xpos = row["start"]
            if needReverse:
                ha = "left"
                genename = "  " + row["name"]
                # xpos = row['end']
            ypos = plot_y + 0.5
            if line == 1:
                xpos = row["start"] + abs(row["start"] - row["end"]) / 2
                ypos = 0.8
                ha = "center"
            ax.text(
                xpos,
                ypos,
                genename + "  ",
                ha=ha,
                va="center",
                color=text_col,
                fontsize=fontszie,
            )
        else:
            ha = "left"
            genename = "  " + row["name"]
            xpos = row["end"]
            if needReverse:
                ha = "right"
                genename = row["name"] + "  "
                # xpos = row['start']

            ypos = plot_y + 0.5
            if line == 1:
                xpos = row["start"] + abs(row["start"] - row["end"]) / 2
                ypos = 0.8
                ha = "center"
            ax.text(
                xpos,
                ypos,
                genename,
                ha=ha,
                va="center",
                color=text_col,
                fontsize=fontszie,
            )

        ii += 1

    xlim_s = start
    xlim_e = end
    if needReverse == True:
        xlim_s = end
        xlim_e = start

    ax.set_xlim(xlim_s, xlim_e)
    ax.set_ylim(top=0, bottom=line)
    if plot_gene_num < line:
        ax.spines["bottom"].set_position(("data", plot_gene_num))

    if ax_on == False:
        spines = ["top", "bottom", "left", "right"]
        for i in spines:
            ax.spines[i].set_visible(False)
            # for i in ['left','top','right']:
            ax.spines[i].set_color("none")
            ax.spines[i].set_linewidth(0)
    ax.spines["bottom"].set_color("black")
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(bottom=True, top=False, left=False, right=False)
    ax.set_xticklabels("")
    ax.set_yticklabels("")
