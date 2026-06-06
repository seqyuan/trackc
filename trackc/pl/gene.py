from typing import Optional, Sequence, Union

import pandas as pd
from matplotlib.axes import Axes

from trackc.pl.bigwig import _make_multi_region_ax
from trackc.tl._getRegionsCmat import GenomeRegion


def gene_track(
    ax: Optional[Axes] = None,
    bed12: Union[pd.DataFrame, str] = None,
    regions: Union[Sequence[str], str, None] = None,
    track_type: Union[str, None] = "gene",
    show_label: Union[bool, Sequence[str], str] = True,
    pos_strand_gene_color: Union[str, None] = "#3366CC",
    neg_strand_gene_color: Union[str, None] = "#EECFA1",
    line: Union[int, None] = 1,
    gene_fontsize: Union[int, None] = 7,
    max_labels: Optional[int] = 60,
    all_labels_inside: bool = False,
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
        Bed12 files can be converted from GTF using `trackc gtf2bed`.
        https://trackc.readthedocs.io/en/latest/analysis_guide/genebed12/bed12.html
    regions: `str` | `str list`
        genome regions, format: `chrom:start-end`.
        e.g. ``['chr18:47950000-48280000', 'chr18:75280000-74850000']`` or ``"chr18:45000000-78077248"``.
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
    gene_fontsize: `int`
        gene label fontsize
    max_labels: `int` | `None`
        Maximum number of genes for automatic label display. The effective
        automatic limit also considers `line` to avoid dense-label plots.
        Set to `None` to disable this limit. Explicit `show_label` gene names
        are not limited.
    all_labels_inside: `bool`
        If True, place labels inside the plotting region when possible.
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
    chrom_names = bed12["chrom"].unique()

    for ix, row in line_GenomeRegions.iterrows():
        if row["chrom"] not in chrom_names:
            raw_chr = row["chrom"]
            if row["chrom"].startswith("chr"):
                row["chrom"] = row["chrom"].lstrip("chr")
            else:
                row["chrom"] = "chr" + row["chrom"]
            if row["chrom"] not in chrom_names:
                print(f"{raw_chr} not in bigwig chroms!")
                return

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
                fontsize=gene_fontsize,
                max_labels=max_labels,
                all_labels_inside=all_labels_inside,
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
    fontsize=5,
    max_labels=60,
    all_labels_inside=False,
    ax_on=False,
):
    gene_bed = gene_bed[gene_bed["chrom"] == chrom]
    region_start = min(start, end)
    region_end = max(start, end)
    gene_bed_plot = gene_bed[
        (gene_bed["start"] <= region_end) & (gene_bed["end"] >= region_start)
    ]
    gene_bed_plot = gene_bed_plot.sort_values(by=["start", "end"])
    # print(gene_bed_plot

    plot_gene_num = gene_bed_plot.shape[0]

    line = max(int(line), 1)
    region_width = abs(end - start)
    head_length = _get_arrow_head_length(region_width, line)
    row_last_position = []
    label_padding = region_width * 0.01
    label_unit = _estimate_label_unit(ax, start, end, fontsize)
    explicit_labels = not isinstance(show_label, bool)
    auto_label_limit = _get_auto_label_limit(max_labels, line)
    display_auto_labels = (
        auto_label_limit is None or explicit_labels or plot_gene_num <= auto_label_limit
    )

    for i, row in gene_bed_plot.iterrows():
        # col = pos_strand_gene_color
        text_col = pos_strand_gene_color

        if row["strand"] == "-":
            # col = neg_strand_gene_color
            text_col = neg_strand_gene_color

        gene_start = max(row["start"], region_start)
        gene_end = min(row["end"], region_end)
        body_left = min(gene_start, gene_end)
        body_right = max(gene_start, gene_end)
        wants_label = _should_show_gene_label(
            show_label, row["name"], display_auto_labels
        )
        label_side, label_left, label_right = _get_gene_label_extent(
            row["name"],
            gene_start,
            gene_end,
            start,
            end,
            needReverse,
            wants_label,
            label_unit,
            label_padding,
            all_labels_inside,
            center_label=line == 1,
        )
        show_this_label = wants_label
        if wants_label:
            layout_left = min(body_left, label_left)
            layout_right = max(body_right, label_right)
            _, label_fits = _get_free_gene_row(
                row_last_position,
                layout_left,
                layout_right,
                line,
                label_padding,
                commit=False,
            )
            if label_fits or explicit_labels:
                plot_y, _ = _get_free_gene_row(
                    row_last_position,
                    layout_left,
                    layout_right,
                    line,
                    label_padding,
                )
            else:
                plot_y, _ = _get_free_gene_row(
                    row_last_position,
                    body_left,
                    body_right,
                    line,
                    label_padding,
                )
                show_this_label = False
        else:
            plot_y, _ = _get_free_gene_row(
                row_last_position,
                body_left,
                body_right,
                line,
                label_padding,
            )

        ax.plot(
            (gene_start, gene_end),
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

        arrow_s = gene_end
        dx = 0.3
        if row["strand"] == "-":
            arrow_s = gene_start
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
        if not show_this_label:
            continue

        label_text = row["name"]
        if label_side == "center":
            ha = "center"
            xpos = gene_start + abs(gene_start - gene_end) / 2
            if line == 1:
                plot_y = 0.3
        elif label_side == "left":
            ha = "right"
            xpos = gene_start - label_padding
            label_text = row["name"] + "  "
        else:
            ha = "left"
            xpos = gene_end + label_padding
            label_text = "  " + row["name"]
        ax.text(
            xpos,
            plot_y + 0.5,
            label_text,
            ha=ha,
            va="center",
            color=text_col,
            fontsize=fontsize,
            clip_on=all_labels_inside,
        )

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


def _estimate_label_unit(ax, start, end, fontsize):
    fig_width = ax.figure.get_figwidth() or 1
    axis_width = ax.get_position().width or 1
    region_width = abs(end - start)
    # Approximate one character as 0.6 em and convert inches to genomic units.
    char_width_in = fontsize / 72 * 0.6
    return region_width * char_width_in / (fig_width * axis_width)


def _get_arrow_head_length(region_width, line):
    if region_width <= 0:
        return 0
    nominal = region_width / max(line + 2, 3) / 8
    return min(nominal, region_width * 0.01)


def _should_show_gene_label(show_label, gene_name, display_auto_labels):
    if isinstance(show_label, bool):
        return show_label and display_auto_labels
    if isinstance(show_label, str):
        return gene_name == show_label
    if isinstance(show_label, (list, tuple, set)):
        return gene_name in show_label
    return bool(show_label)


def _get_gene_label_extent(
    gene_name,
    gene_start,
    gene_end,
    region_start,
    region_end,
    need_reverse,
    show_label,
    label_unit,
    label_padding,
    all_labels_inside,
    center_label=False,
):
    label_width = (len(str(gene_name)) + 2) * label_unit if show_label else 0
    region_left = min(region_start, region_end)
    region_right = max(region_start, region_end)
    gene_left = min(gene_start, gene_end)
    gene_right = max(gene_start, gene_end)

    if not show_label:
        return "none", gene_left, gene_right + 2 * label_padding

    if center_label or (
        label_width >= max(gene_right - gene_left, label_padding) and all_labels_inside
    ):
        center = gene_left + (gene_right - gene_left) / 2
        return "center", center - label_width / 2, center + label_width / 2

    label_side = "left" if need_reverse else "right"
    if label_side == "right":
        label_left = gene_right
        label_right = gene_right + label_padding + label_width
        if all_labels_inside and label_right > region_right:
            left_start = gene_left - label_padding - label_width
            if left_start >= region_left:
                return "left", left_start, gene_right
            center = gene_left + (gene_right - gene_left) / 2
            return "center", center - label_width / 2, center + label_width / 2
    else:
        label_left = gene_left - label_padding - label_width
        label_right = gene_left
        if all_labels_inside and label_left < region_left:
            right_end = gene_right + label_padding + label_width
            if right_end <= region_right:
                return "right", gene_left, right_end
            center = gene_left + (gene_right - gene_left) / 2
            return "center", center - label_width / 2, center + label_width / 2

    return label_side, label_left, label_right


def _get_auto_label_limit(max_labels, line):
    if max_labels is None:
        return None
    return min(int(max_labels), max(int(line) * 4, 4))


def _get_free_gene_row(row_last_position, left, right, max_rows, padding, commit=True):
    for row_idx, row_end in enumerate(row_last_position):
        if left > row_end + padding:
            if commit:
                row_last_position[row_idx] = right
            return row_idx, True
    if len(row_last_position) < max_rows:
        if commit:
            row_last_position.append(right)
        return len(row_last_position), True
    row_idx = min(range(max_rows), key=lambda idx: row_last_position[idx])
    if commit:
        row_last_position[row_idx] = max(row_last_position[row_idx], right)
    return row_idx, False
