from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
import pandas as pd
import numpy as np
from trackc.tl._getRegionsCmat import GenomeRegion
from .bigwig import _make_multi_region_ax


def _two_degree_bc(x_l=10, x_r=90, y_lr=0, y2=10, dots_num=100): 
    """
        bezier curve for loop links
    """

    xt = []
    yt = []
    x_mid = (x_l + x_r)/2
    x_dots12 = np.linspace(x_l, x_mid, dots_num)
    y_dots12 = np.linspace(y_lr, y2, dots_num)
    x_dots23 = np.linspace(x_mid, x_r, dots_num)
    y_dots23 = np.linspace(y2, y_lr, dots_num)
    for i in range(dots_num):
        x = x_dots12[i] + (x_dots23[i]-x_dots12[i])*i / (dots_num-1)
        y = y_dots12[i] + (y_dots23[i]-y_dots12[i])*i / (dots_num-1)
        xt.append(x)
        yt.append(y)
    return (xt, yt)

def _plot_loop_arc(ax, loop_df, color, max_extend, invert_y, start, end, left_anchor, right_anchor):
    if loop_df.shape[0] == 0:
        return
    
    top_y = 0

    for i, row in loop_df.iterrows():
        top = row['length']/max_extend
        if top < 0.5:
            top = 0.5
        elif top >0.8:
            top = 0.8
        else:
            pass
        
        xt, yt = _two_degree_bc(x_l=row[left_anchor], x_r=row[right_anchor], y_lr=0, y2=top, dots_num=100)
    
        ax.plot(xt, yt, color=color, linewidth=0.5, solid_capstyle='butt')
        if max(yt) > top_y:
            top_y = max(yt)
            
    ax.set_xlim(start, end)
    if invert_y==True:
        ax.set_ylim(0.5, 0)
    else:
        ax.set_ylim(0, 0.5)



def links_track(
        ax: Optional[Axes] = None,
        data: pd.DataFrame = None, 
        regions: Union[Sequence[str], str, None] = None,
        links_type: Union[str, None] = 'arc',
        color: Union[Sequence[str], None] = '#66AC84',
        maxrange: Union[int, None] = 3000000,
        invert_y: Optional[bool] = False,
        anchor: Union[str, None] = 'inside',
        label: Optional[str] = None,
        label_rotation: Union[int, None] = 0,
        label_fontsize: Optional[int] = 12,
        ax_on: bool = False,
        ):
    """\
    Plot loop arc, support for multiple or reverse genome regions.
    
    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    data: `pd.DataFrame`
        the file format expected is like this:
        chr1 x1 x2 chr2 y1 y2
        the fields after the y2 will be ignored,
        recommend the result of juicer loop calling
    regions: `str` | `str list`
        The genome regions to show the arc.
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means the reverse region
    links_type: `str`
        Optional is ['arc']
    color: `str` or `str list`
        the color of the arc links, if multiple regions, color para can  set as list
    maxrange: `int`
        The maximum distance between two anchor of loop, to filter loops
        if value is None, all loop will bed plotted 
    invert_y: `bool`
        whether reverse the y-axis
    anchor: `str`
        Optional is ['inside', 'outside']
        inside: arc link x2 y1
        outside: arc link x2 y1
    label: `str`
        the title of the track, will show on the left
    label_rotation: `int`
        the label text rotation
    label_fontsize: `int`
        the label text fontsize
    ax_on: `bool`
        whether show the spines
        
    Example
    -------
    >>> import trackc as tc
    >>> regions = ['chr7:153000000-151000000', 'chr11:118500000-116500000']

    >>> ten = tc.tenon(width=8, height=1)
    >>> ten.add(pos='bottom', height=1)
    >>> ten.add(pos='bottom', height=1, hspace=0.1)
    >>> ten.add(pos='bottom', height=0.4, hspace=0.1)

    >>> tc.pl.links_track(ax=ten.axs(0), data=loops, label='GM12878', regions=regions, 
                color=['#66AC84', 'tab:purple'], maxrange=3000000, anchor='inside')
    >>> tc.pl.links_track(ax=ten.axs(1), data=loops, label='GM12878', regions=regions, 
                color='tab:purple', invert_y=True, anchor='outside', ax_on=True)
    >>> tc.pl.multi_scale_track(ten.axs(2), regions=regions, scale_adjust='Mb', intervals=1, tick_rotation=0, tick_fontsize=10, colors=['#66AC84', 'tab:purple'])
    >>> tc.savefig('trackc_links_track.pdf')
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    axs = _make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    if isinstance(color, list)==False:
        color = [color]
    if len(color) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(color) - 1) // len(color)
        color = (color * repeat_times)[:line_GenomeRegions.shape[0]]
    
    data = data.iloc[:,[0,1,2,3,4,5]]
    data.columns = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

    if anchor == "mid":
        left_anchor = 'left_mid'
        right_anchor = 'right_mid'

    if anchor == "inside":
        left_anchor = 'x2'
        right_anchor = 'y1'
    if anchor == "outside":
        left_anchor = 'x1'
        right_anchor = 'y2'

    max_extend = 0
    links_df_list = []

    data['chr1'] = data['chr1'].astype(str)
    data['chr2'] = data['chr2'].astype(str)

    for ix, row in line_GenomeRegions.iterrows(): 
        loop_bed_plot = data[data['chr1']==row['chrom']]
        

        loop_bed_plot['length'] = loop_bed_plot["y2"] - loop_bed_plot["x1"]
        
        loop_bed_plot = loop_bed_plot[((loop_bed_plot['x1'] >= row['fetch_start']) & (loop_bed_plot['x1'] <= row['fetch_end'])) | \
                                ((loop_bed_plot['x2'] >= row['fetch_start']) & (loop_bed_plot['x2'] <= row['fetch_end']) |\
                                (loop_bed_plot['y1'] >= row['fetch_start']) & (loop_bed_plot['y1'] <= row['fetch_end'])) | \
                                    ((loop_bed_plot['y2'] >= row['fetch_start']) & (loop_bed_plot['y2'] <= row['fetch_end']))]
        
        if isinstance(maxrange, int):
            loop_bed_plot = loop_bed_plot[(loop_bed_plot["y2"] - loop_bed_plot["x1"])<maxrange] 
        #print('looop num:', loop_bed_plot.shape)
        if loop_bed_plot.shape[0] == 0:
            #for i in ['top', 'right', "left", "bottom"]:
            #    ax.spines[i].set_color('none')
            #    ax.spines[i].set_linewidth(0)
            next

        loop_bed_plot['left_mid'] = (loop_bed_plot['x1']+loop_bed_plot['x2'])/2
        loop_bed_plot['right_mid'] = (loop_bed_plot['y1']+loop_bed_plot['y2'])/2

        links_df_list.append(loop_bed_plot)
        if loop_bed_plot.shape[0] >0:
            if max_extend < max(loop_bed_plot['length']):
                max_extend = max(loop_bed_plot['length'])
   
    for ix, row in line_GenomeRegions.iterrows(): 
        if links_type == 'arc':
            _plot_loop_arc(axs[ix], links_df_list[ix], color[ix], max_extend, invert_y, row['start'], row['end'], left_anchor, right_anchor)
        #if links_type == 'loop':
            #plot_triangle()

    ax.set_ylabel(label, fontsize=label_fontsize, rotation=label_rotation, 
                  horizontalalignment='right',verticalalignment='center')
    
    
    if ax_on == False:
        spines = ['top', 'bottom', 'left', 'right']
        if invert_y == True:
            del spines[0]
        else:
            del spines[1]
        for i in spines:
            ax.spines[i].set_visible(False)

    ax.set_xticks([])
    ax.set_xticklabels('')
    ax.set_yticks([])
    ax.set_yticklabels('')

