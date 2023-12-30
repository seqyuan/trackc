import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

from typing import List, Union

from matplotlib.axes import Axes


def make_spec(
    figsize: tuple = (8, 3),
    height_ratios: Union[List[float], None] = None,
    width_ratios: Union[List[float], None] = None,
    hspace: float = 0,
    wspace: float = 0,
):
    """
    Make `GridSpec` allows complex layout of Axes in the figure.

    Parameters
    ----------
    figsize : 2-tuple of floats, default: `figure.figsize`
            Figure dimension ``(width, height)`` in inches.
            ``Whole Figure Size``
    width_ratios : `array-like`
        if `height_ratios` is setting, this parameter will ignore.
        Defines the relative widths of the columns. Each column gets a
        relative width of ``width_ratios[i] / sum(width_ratios)``.
        If not given, all columns will have the same width.
    height_ratios : `array-like`
        Defines the relative heights of the rows. Each row gets a
        relative height of ``height_ratios[i] / sum(height_ratios)``.
        If not given, all rows will have the same height.
    wspace, hspace : `float`, default: 0
        The amount of width/height reserved for space between subfigures,
        expressed as a fraction of the average subfigure width/height.
        If not given, the values will be inferred from a figure or
        rcParams when necessary.
    Returns
    -------
    Figure, list[Axes]
        `Figure`: top level `~matplotlib.artist.Artist`, which holds all plot elements.
            Many methods are implemented in `FigureBase`.
        `list[Axes]`: is a list of the Axes objects.  The order of the axes is left-to-right
            or top-to-bottom of their position in the total layout.

    Example
    -------
    >>> import trackc as tc
    >>> fig, axs = tc.make_spec(figsize=(6,6), height_ratios=[1,1,1])
    >>> axs[0].plot([1,2,3])
    >>> tc.savefig('test.pdf')
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
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        return fig, ax

    fig = plt.figure(figsize=figsize)
    spec = fig.add_gridspec(
        nrows=nrows,
        ncols=ncols,
        left=None,
        bottom=None,
        right=None,
        top=None,
        wspace=wspace,
        hspace=hspace,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )

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
    """
    Make a virtual-figure, based on virtual-figure, users can `add` `Axes` to
    the `top or bottom` of the total layout. After adding `Axes`, the user can use the ``axs`` method
    to select an axes from the added axes list. The order of the axes is top-to-bottom of
    the position in the total layout.

    Just like Lego blocks, you can add new blocks above (top) or below (bottom)

    Example
    -------
    >>> import trackc as tc
    >>> ten = trackc.tenon(figsize=(6,1))
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='top', height=1, hspace=0.1)
    >>> ten.axs(0).plot([1,2,3])
    >>> ten.axs(1).plot([1,1])
    >>> tc.savefig('test.pdf')
    """

    mortises = []
    fig = None
    ax = None

    def __init__(self, figsize=(7, 1)):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_axis_off()
        self.fig = fig
        self.ax = ax

    def show(self):
        for i, v in enumerate(self.mortises):
            # v.ax.set_ylabel('mortises[{0}]'.format(i), rotation=0)
            v.ax.text(0.5, 0.5, ".axs({0})".format(i), ha="center", va="center")
            v.ax.set_xticklabels([])
            v.ax.set_xticks([])
        self.fig.show()

    def add(self, pos="top", height=1, hspace=0):
        """
        Add a `Axes` to the total layout.
        Just like Lego blocks, you can add new blocks above (top) or below (bottom)

        Parameters
        ----------
        pos
            one option of ['top', 'bottom']
            Choose whether to add a new Axes at the top or bottom.
            If set top, then a `~matplotlib.axes.Axes` object will add to top of the top `~matplotlib.axes.Axes`
        height: `float`
            the relative height of the newly added `Axes` compared to the virtual-figure.
            For example, if the virtual-figure has a height of 1 and the new Axes has a height of 4,
            the actual height of the subplot will be 1 * 4.
        hspace: `float`
            The amount of height reserved for space between the top or bottom subfigures
        """
        if len(self.mortises) != 0:
            ix = -1
            if pos == "top":
                ix = 0
            ref_lego_ax = self.mortises[ix]
        else:
            ref_lego_ax = mortise(self.ax, 0, height, 0)
            self.mortises = [ref_lego_ax]
            return

        bottom = ref_lego_ax.bottom - hspace - height
        if pos == "top":
            bottom = ref_lego_ax.top + hspace

        lego_ax_demo = mortise(self.ax, bottom, height, hspace)
        if pos == "top":
            self.mortises.insert(0, lego_ax_demo)
        else:
            self.mortises.append(lego_ax_demo)

    def _remove(self, pos="top"):
        """
        Remove a `Axes` to the total layout.
        Just like Lego blocks, you can remove new blocks from above (top) or below (bottom)

        Parameters
        ----------
        pos
            Option of ['top', 'bottom']
            Choose whether to remove a Axes from top or bottom of the whole layout.
        """
        if len(self.mortises) >= 1:
            if pos == "top":
                self.mortises.remove(self.mortises[0])
            elif pos == "bottom":
                self.mortises.remove(self.mortises[-1])
            else:
                pass

    def axs(self, index):
        """
        Get one Axes from the total layout.

        Parameters
        ----------
        index
            The Axes list element index

        Returns
        -------
        `~matplotlib.axes.Axes` object
        """
        return self.mortises[index].ax


def savefig(outfile):
    """
    Save the figure.
    Parameters
    ----------
    outfile
        out file path
    """
    # plt.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")
    # plt.savefig(outfile)


def show():
    """
    show the Axes of the layout, help users resize each axes.
    """
    plt.tight_layout()
