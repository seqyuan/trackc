import ast
import sys
import warnings
from typing import Any, Callable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyBigWig
from matplotlib.axes import Axes

from trackc.tl._getRegionsCmat import GenomeRegion

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
    number_of_bins: Optional[int] = None,
    style: Optional[str] = "bar",
    summary_type: Union[str, None] = "mean",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    primary_col: Union[Sequence[str], None] = "#3271B2",
    secondary_col: Union[Sequence[str], None] = "#FBD23C",
    alpha: Optional[float] = 1.0,
    nans_to_zeros: bool = True,
    transform: str = "no",
    log_pseudocount: Optional[float] = 0.0,
    invert_y: Optional[bool] = False,
    label: Optional[str] = None,
    label_rotation: Union[int, None] = 0,
    label_fontsize: Optional[int] = 9,
    tick_fontsize: Optional[int] = 7,
    tick_fl: Optional[str] = "%0.2f",
    ax_on: bool = False,
    second_bw=None,
    operation: str = "file",
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
        binsize divided to computing signal summary statistics.
        Ignored when ``number_of_bins`` is provided.
    number_of_bins: `int`
        Number of bins used to summarize each region.
    style: `str`
        plot type, default='bar', options=['bar', 'line', 'fill']
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
    nans_to_zeros: `bool`
        If True, missing values are converted to zero before plotting.
    transform: `str`
        Signal transform, one of ['no', 'log', 'log1p', '-log', 'log2', 'log10']
    log_pseudocount: `float`
        Pseudocount used by log-based transforms.
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
    second_bw: `pyBigWig.open` query object, or bigwig file path
        Optional second bigWig used by ``operation``.
    operation: `str`
        Operation applied to the summarized values. The first bigWig is exposed
        as ``file`` and the optional second bigWig as ``second_file``.
        Examples: ``0.89 * file``, ``file - second_file``,
        ``log2((1 + file) / (1 + second_file))``.
        ``operation`` and ``transform`` cannot be used together.

    Example
    -------
    >>> import trackc as tc
    >>> import pyBigWig
    >>> H3K27ac = pyBigWig.open('./GSM4604189.bigwig')

    >>> ten = tc.tenon(figsize=(8,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=1, hspace=0.2)

    >>> regions = ['chr8:127000000-129200000', 'chr14:96500000-99300000']
    >>> tc.pl.bw_track(H3K27ac, ten.axs(0), regions=regions, vmax=20, label='H3K27ac', binsize=10000, primary_col='tab:blue')
    >>> tc.pl.bw_track(H3K27ac, ten.axs(1), regions=regions, vmax=5, label='H3K27ac', binsize=10000, invert_y=True, ax_on=True)
    >>> tc.savefig('trackc_bigwig_track.pdf')
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat(
            [GenomeRegion(i).GenomeRegion2df() for i in regions]
        )
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    bw = _open_bigwig_if_needed(bw)
    second_bw = _open_bigwig_if_needed(second_bw)
    operation, compiled_operation, needs_second_bw = _compile_operation(operation)

    if operation != "file" and transform != "no":
        raise ValueError("'operation' and 'transform' cannot be set at the same time")
    if needs_second_bw and second_bw is None:
        raise ValueError(
            f"operation '{operation}' requires a second bigWig via second_bw"
        )

    axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    if isinstance(primary_col, list) == False:
        primary_col = [primary_col]
    if len(primary_col) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(primary_col) - 1) // len(
            primary_col
        )
        primary_col = (primary_col * repeat_times)[: line_GenomeRegions.shape[0]]

    if isinstance(secondary_col, list) == False:
        secondary_col = [secondary_col]
    if len(secondary_col) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(secondary_col) - 1) // len(
            secondary_col
        )
        secondary_col = (secondary_col * repeat_times)[: line_GenomeRegions.shape[0]]

    min_y = 0
    max_y = 0
    plot_bottom_line = True
    has_finite_values = False

    for i, row in line_GenomeRegions.iterrows():
        bins = _resolve_nbins(row["len"], binsize, number_of_bins)
        plot_values = _get_bigwig_stats_array(
            bw,
            row["chrom"],
            row["fetch_start"],
            row["fetch_end"],
            bins,
            summary_type,
        )
        if plot_values is None:
            return

        second_values = None
        if needs_second_bw:
            second_values = _get_bigwig_stats_array(
                second_bw,
                row["chrom"],
                row["fetch_start"],
                row["fetch_end"],
                bins,
                summary_type,
            )
            if second_values is None:
                return

        if nans_to_zeros:
            plot_values = np.nan_to_num(plot_values, nan=0.0)
            if second_values is not None:
                second_values = np.nan_to_num(second_values, nan=0.0)

        plot_values = _apply_operation(
            plot_values, second_values, operation, compiled_operation
        )
        plot_values = _transform_scores(plot_values, transform, log_pseudocount)
        if np.any(np.isfinite(plot_values) & (plot_values < 0)):
            plot_bottom_line = False

        x_values = np.arange(plot_values.shape[0], dtype=float)
        _plot_bigwig_values(
            axs[i],
            x_values,
            plot_values,
            style,
            primary_col[i],
            secondary_col[i],
            alpha,
        )

        right, left = plot_values.shape[0], 0
        if row["isReverse"] == True:
            left, right = bins, 0
        axs[i].set_xlim(left, right)

        finite_values = plot_values[np.isfinite(plot_values)]
        if finite_values.size > 0:
            has_finite_values = True
            if min_y > float(np.min(finite_values)):
                min_y = float(np.min(finite_values))
            if max_y < float(np.max(finite_values)):
                max_y = float(np.max(finite_values))

    if vmin == None:
        vmin = min_y if has_finite_values else 0
    if vmax == None:
        vmax = max_y if has_finite_values else 0

    if vmin == vmax:
        if vmin == 0:
            vmax = 1
        else:
            delta = abs(vmin) * 0.05
            vmin -= delta
            vmax += delta

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
        ax.set_yticklabels(
            [f"{tick_fl % vmin}", "0", f"{tick_fl % vmax}"], fontsize=tick_fontsize
        )
        if invert_y:
            ax.set_yticks([vmax, 0, vmin])
            ax.set_yticklabels(
                [f"{tick_fl % vmax}", "0", f"{tick_fl % vmin}"], fontsize=tick_fontsize
            )

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
            if plot_bottom_line == True:
                del spines[0]
            else:
                del spines[2]
        else:
            if plot_bottom_line == True:
                del spines[1]
            else:
                del spines[2]
        for i in spines:
            ax.spines[i].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels("")


def _resolve_nbins(region_len, binsize, number_of_bins):
    if number_of_bins is not None:
        bins = int(number_of_bins)
        if bins < 1:
            raise ValueError("number_of_bins must be >= 1")
        return bins

    if binsize is None:
        raise ValueError("binsize or number_of_bins must be provided")

    bins = int(region_len / binsize)
    return max(bins, 1)


def _open_bigwig_if_needed(bw):
    if isinstance(bw, str):
        return pyBigWig.open(bw)
    return bw


def _get_bigwig_stats_array(bw, chrom, fetch_start, fetch_end, bins, summary_type):
    query_chrom = _resolve_bigwig_chrom_name(bw, chrom)
    if query_chrom is None:
        print(f"{chrom} not in bigwig chroms!")
        return None

    plot_list = bw.stats(
        query_chrom,
        int(fetch_start),
        int(fetch_end),
        type=summary_type,
        nBins=bins,
    )
    return np.array([np.nan if v is None else float(v) for v in plot_list], dtype=float)


def _resolve_bigwig_chrom_name(bw, chrom):
    if chrom in bw.chroms():
        return chrom

    if chrom.startswith("chr"):
        alt_chrom = chrom.lstrip("chr")
    else:
        alt_chrom = "chr" + chrom

    if alt_chrom in bw.chroms():
        return alt_chrom
    return None


def _compile_operation(operation):
    if operation is None:
        operation = "file"

    operation = str(operation).strip()
    if operation == "":
        raise ValueError("operation must not be empty")
    if operation == "file":
        return operation, None, False

    try:
        tree = ast.parse(operation, mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"invalid operation: {operation}") from exc

    allowed_functions = {
        "abs",
        "clip",
        "exp",
        "log",
        "log1p",
        "log2",
        "log10",
        "max",
        "min",
        "sqrt",
        "where",
    }
    allowed_names = {"file", "second_file"} | allowed_functions
    allowed_nodes = (
        ast.Expression,
        ast.BinOp,
        ast.UnaryOp,
        ast.Call,
        ast.Name,
        ast.Load,
        ast.Constant,
        ast.Add,
        ast.Sub,
        ast.Mult,
        ast.Div,
        ast.Pow,
        ast.Mod,
        ast.UAdd,
        ast.USub,
    )

    for node in ast.walk(tree):
        if not isinstance(node, allowed_nodes):
            raise ValueError(
                "operation contains unsupported syntax; use arithmetic with "
                "file, second_file, and supported functions"
            )
        if isinstance(node, ast.Name) and node.id not in allowed_names:
            raise ValueError(f"unsupported name in operation: {node.id}")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in allowed_functions:
                raise ValueError("operation contains an unsupported function call")
            if len(node.keywords) > 0:
                raise ValueError("operation does not support keyword arguments")

    needs_second_bw = any(
        isinstance(node, ast.Name) and node.id == "second_file"
        for node in ast.walk(tree)
    )
    return operation, compile(tree, "<bw_track operation>", "eval"), needs_second_bw


def _apply_operation(file_values, second_values, operation, compiled_operation):
    file_values = np.array(file_values, dtype=float, copy=True)
    if operation == "file":
        return file_values

    if second_values is None and "second_file" in operation:
        raise ValueError(
            f"operation '{operation}' requires a second bigWig via second_bw"
        )

    namespace = {
        "abs": np.abs,
        "clip": np.clip,
        "exp": np.exp,
        "file": file_values,
        "log": np.log,
        "log1p": np.log1p,
        "log2": np.log2,
        "log10": np.log10,
        "max": _operation_max,
        "min": _operation_min,
        "sqrt": np.sqrt,
        "where": np.where,
    }
    if second_values is not None:
        namespace["second_file"] = np.array(second_values, dtype=float, copy=True)

    try:
        result = eval(compiled_operation, {"__builtins__": {}}, namespace)
    except Exception as exc:
        raise ValueError(f"failed to evaluate operation '{operation}': {exc}") from exc

    result = np.asarray(result, dtype=float)
    if result.shape == ():
        result = np.full(file_values.shape, float(result), dtype=float)
    if result.shape != file_values.shape:
        raise ValueError(
            "operation must return one value per bin, matching the first bigWig"
        )
    return result


def _operation_max(*values):
    if len(values) == 0:
        raise ValueError("max requires at least one argument")

    result = np.asarray(values[0], dtype=float)
    for value in values[1:]:
        result = np.maximum(result, np.asarray(value, dtype=float))
    return result


def _operation_min(*values):
    if len(values) == 0:
        raise ValueError("min requires at least one argument")

    result = np.asarray(values[0], dtype=float)
    for value in values[1:]:
        result = np.minimum(result, np.asarray(value, dtype=float))
    return result


def _transform_scores(scores, transform, log_pseudocount):
    transformed = np.array(scores, dtype=float, copy=True)
    finite_mask = np.isfinite(transformed)
    if not finite_mask.any() or transform == "no":
        return transformed

    min_value = float(np.nanmin(transformed))
    if transform in ["log", "log2", "log10"]:
        if min_value <= -log_pseudocount:
            raise ValueError(
                f"{transform} transform requires all values > {-log_pseudocount}"
            )
        ops = {"log": np.log, "log2": np.log2, "log10": np.log10}
        transformed[finite_mask] = ops[transform](
            log_pseudocount + transformed[finite_mask]
        )
        return transformed

    if transform == "log1p":
        if min_value <= -1:
            raise ValueError("log1p transform requires all values > -1")
        transformed[finite_mask] = np.log1p(transformed[finite_mask])
        return transformed

    if transform == "-log":
        if min_value <= -log_pseudocount:
            raise ValueError(f"-log transform requires all values > {-log_pseudocount}")
        transformed[finite_mask] = -np.log(log_pseudocount + transformed[finite_mask])
        return transformed

    raise ValueError(
        "transform must be one of ['no', 'log', 'log1p', '-log', 'log2', 'log10']"
    )


def _plot_bigwig_values(
    ax,
    x_values,
    plot_values,
    style,
    primary_col,
    secondary_col,
    alpha,
):
    finite_mask = np.isfinite(plot_values)
    pos_mask = finite_mask & (plot_values >= 0)
    neg_mask = finite_mask & (plot_values < 0)

    if style == "line":
        if primary_col == secondary_col:
            ax.plot(x_values, plot_values, color=primary_col, alpha=alpha)
        else:
            pos_values = np.where(pos_mask, plot_values, np.nan)
            neg_values = np.where(neg_mask, plot_values, np.nan)
            ax.plot(x_values, pos_values, color=primary_col, alpha=alpha)
            ax.plot(x_values, neg_values, color=secondary_col, alpha=alpha)
        return

    if style == "fill":
        pos_values = np.where(pos_mask, plot_values, np.nan)
        neg_values = np.where(neg_mask, plot_values, np.nan)
        ax.fill_between(
            x_values,
            0,
            pos_values,
            facecolor=primary_col,
            color=primary_col,
            linewidth=0,
            alpha=alpha,
        )
        ax.fill_between(
            x_values,
            0,
            neg_values,
            facecolor=secondary_col,
            color=secondary_col,
            linewidth=0,
            alpha=alpha,
        )
        return

    ax.bar(
        x=x_values[pos_mask],
        height=plot_values[pos_mask],
        width=1,
        bottom=0,
        color=primary_col,
        align="edge",
        edgecolor=None,
        alpha=alpha,
    )
    ax.bar(
        x=x_values[neg_mask],
        height=plot_values[neg_mask],
        width=1,
        bottom=0,
        color=secondary_col,
        align="edge",
        edgecolor=None,
        alpha=alpha,
    )
