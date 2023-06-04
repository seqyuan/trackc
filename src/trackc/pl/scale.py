from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
import pandas as pd
from trackc.tl import GenomeRegion
from trackc.pa import trackcl_11

def _scale_ticks(ticks, scale='bp', tick_fl='%0.2f'):
    if scale == 'Mb':
        ticks = ticks/1000000
        ticks = ["{0}".format(tick_fl % i) for i in ticks]
    elif scale == 'kb': 
        ticks = ticks/1000
        ticks = ["{0}".format(tick_fl % i) for i in ticks]
    else:
        return ticks
    return ticks
def scale_track(ax: Optional[Axes] = None,
                region: Union[str, None] = None,
                tick_pos: str = 'bottom',
                ratio2ax: float = 0.5,
                label_fontsize: Union[int, None] = 10,
                scale_adjust: Union[str, None] = 'kb',
                tick_fl: Union[str, None] ='%0.2f',
                tick_fontsize: Union[int, None] = 8,
                tick_rotation: Union[int, None] = 0,
                space: float = 0.1,
                ):
    """\
    Plot one region scale bar track

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    region: `str`
        one genome region, format: `chrom:start-end`
        examples: "chr18:45000000-78077248"
        if the start is bigger than end, the genome region will be reversed
    tick_pos: `str`
        ticks position, can be one of ['top', 'bottom']
    ratio2ax: `float`
        the height ratio refer to the given `ax`'s height
    label_fontsize: `int`
        the region text fontsize
    scale_adjust: `str`
        adjust the scale unit to make it pretty, can be one of ['kb', 'Mb']   
    tick_fl: `str`  
        ticks retains a few decimal places
    tick_fontsize: `int`
        ticks text fontsize
    tick_rotation: `int`
        ticks text rotation
    space: `float`
        space relative to the ax
    
    Example
    -------
    >>> import trackc as tc
    >>> region = 'chr7:153000000-151000000'
    
    >>> ten = tc.tenon(width=8, height=1.2)
    >>> ten.add(pos='bottom', height=0.5)
    >>> ten.axs(0).axis('off')

    >>> tc.pl.scale_track(ax=ten.axs(0), region=region, scale_adjust='Mb', tick_pos='bottom', ratio2ax=1.2)
    >>> tc.pl.scale_track(ax=ten.axs(0), region=region, scale_adjust='Mb', tick_pos='top', ratio2ax=1.2)
    >>> tc.savefig('trackc_scalebar_track.pdf')
    """

    line_GenomeRegions = None
    if isinstance(region, list):
        print('scale_track is only for one region')
        return
    else:
        line_GenomeRegions = GenomeRegion(region).GenomeRegion2df()

    line_GenomeRegions['raw_region'] = line_GenomeRegions.index
    line_GenomeRegions = line_GenomeRegions.reset_index()

    chrom = line_GenomeRegions.loc[0, "chrom"]
    start = line_GenomeRegions.loc[0, "start"]
    end = line_GenomeRegions.loc[0, "end"]
    raw_region = line_GenomeRegions.loc[0, 'raw_region']

    pos_dic = {
        'left': [-ratio2ax,	0, ratio2ax, 1],
        'right': [1, 0,	ratio2ax, 1],
        'top': [0, 1+space, 1, ratio2ax],
        'bottom': [0, -ratio2ax-space, 1, ratio2ax]
    }

    ax2 = ax.inset_axes(pos_dic[tick_pos], facecolor='none')
    ax = ax2
    
    if tick_pos in ['top', 'bottom']:
        ax.set_xlim([start, end])

        xticks = ax.get_xticks()
        xtick_labels = xticks

        if scale_adjust == 'Mb':
            xtick_labels = xtick_labels/1000000
            xtick_labels = ["{0}".format(tick_fl % i) for i in xtick_labels]
            ax.set_xticks(xticks, xtick_labels, fontsize=tick_fontsize, rotation=tick_rotation)
            ax.spines['bottom'].set_position(('data', 0))
            #ax.text(end-(abs(end-start)*0.05), -1, 'Mb', fontsize=chrom_fontsize)
        elif scale_adjust == 'kb':
            xtick_labels = xtick_labels/1000
            xtick_labels = ["{0}".format(tick_fl % i) for i in xtick_labels]
            ax.set_xticks(xticks, xtick_labels)
            #ax.text(end-(abs(end-start)*0.05), -1, 'kb', fontsize=chrom_fontsize)
            ax.spines['bottom'].set_position(('data', 0))
        else:
            pass

        spines = ['bottom', 'top', 'right', 'left']
        #chrom_x = start + (end-start)/2
        chrom_x = start
        chrom_y = 1
        va = 'top'
        ha2  = 'right'
        if tick_pos=='bottom':
            del spines[1]
            ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
            chrom_y = 0
            va = 'bottom'
            
        if tick_pos=='top':
            ha2  = 'left'
            del spines[0]
            ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
        for i in spines:
            ax.spines[i].set_visible(False)

        ax.text(chrom_x, chrom_y, raw_region, fontsize=label_fontsize, ha='left', va=va)
        labels = [label.get_text() for label in ax.get_xticklabels()]
        labels[-1] += '({0})'.format(scale_adjust)
        ax.set_xticklabels(labels, fontsize=tick_fontsize, rotation=tick_rotation, va='center', ha=ha2)

        ax.tick_params(which='major', direction='in', pad=-16) 
        ax.set_yticks([])
        ax.set_yticklabels('')  
    


def multi_scale_track(ax: Optional[Axes] = None,
                      regions: Union[Sequence[str], str, None] = None,
                      colors: Union[Sequence[str], None] = None,
                      alpha: Union[float, None] = 1,
                      intervals: Union[int, None] = 1,
                      scale_adjust: Union[str, None] = 'kb',
                      tick_fl: Union[str, None] ='%0.2f',
                      tick_fontsize: Union[int, None] = 8,
                      tick_rotation: Union[int, None] = 0,
                ):
    """\
    Plot region scale bar, support for multiple or reverse genome regions.

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    regions: `str list` | `str`
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means you want to get the reverse region
    colors: `str list`
        scale bar colors
    alpha: `float`
        scale bar alpha
    intervals: `int`
        rows of the scale bar by region
    scale_adjust: `str`
        options in ['kb', 'Mb']
    tick_fl: `str`  
        ticks retains a few decimal places
    tick_fontsize: `int`
        ticks text fontsize
    tick_rotation: `int`
        ticks text rotation

    Example
    -------
    >>> import trackc as tc
    >>> regions = ['7:153000000-151000000', '11:118500000-116500000']
    >>> ten = tc.tenon(width=8, height=1)
    >>> ten.add(pos='bottom', height=1)
    >>> tc.pl.multi_scale_track(ten.axs(0), regions=regions, scale_adjust='Mb', intervals=2)
    """

    if colors != None:
        track_colors  = colors
    else:
        track_colors = trackcl_11

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    line_GenomeRegions = line_GenomeRegions.reset_index()

    if len(track_colors) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(track_colors) - 1) // len(track_colors)
        track_colors = (track_colors * repeat_times)[:line_GenomeRegions.shape[0]]
    else:
        track_colors = track_colors[:line_GenomeRegions.shape[0]]
    line_GenomeRegions['color'] = track_colors
    line_GenomeRegions['len'] = line_GenomeRegions['fetch_end'] - line_GenomeRegions['fetch_start']
    line_GenomeRegions['acum'] = line_GenomeRegions['len'].cumsum()
    sum_len = line_GenomeRegions['len'].sum()
    bottom = list(range(intervals))  
    # broadcast
    if intervals < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + intervals - 1) // intervals
        bottom = (bottom * repeat_times)[:line_GenomeRegions.shape[0]]
    else:
        bottom = bottom[:line_GenomeRegions.shape[0]]
    line_GenomeRegions.loc[:, 'bottom'] = bottom

    line_GenomeRegions['tick_s'] = _scale_ticks(line_GenomeRegions['start'], scale=scale_adjust, tick_fl=tick_fl)
    line_GenomeRegions['tick_e'] = _scale_ticks(line_GenomeRegions['end'], scale=scale_adjust, tick_fl=tick_fl)
    if scale_adjust in ['Mb', 'kb']:
        if intervals > 1:
            line_GenomeRegions.loc[line_GenomeRegions.index[-1], "tick_e"] = line_GenomeRegions.loc[line_GenomeRegions.index[-1], "tick_e"] + "(" + scale_adjust + ")"
        else:
            line_GenomeRegions['tick_e'] = line_GenomeRegions['tick_e'] + "(" + scale_adjust + ")"

    for i, row in line_GenomeRegions.iterrows():
        arrow_s = row['acum'] - row['len']
        dx = row['len']
        if row["isReverse"] == True:
            arrow_s = row['acum']
            dx = -1 * row['len']
            
        ax.arrow(arrow_s, row['bottom'], dx, 0, 
            overhang=1, width=0.02,
            head_width=0.3,
            head_length=sum_len/100,
            length_includes_head=True,
            color=row['color'],
            linewidth=1,
            alpha = alpha
            )
        
        tick_ha = 'center'
        if tick_rotation != 0:
            tick_ha = 'right'
        if intervals > 1:
            ax.text(row['acum'] - row['len']/2, row['bottom']+0.1, row['chrom'], ha='center', va='bottom', fontsize=tick_fontsize, color=row['color'])
            ax.text(row['acum'] - row['len'], row['bottom']-0.1, row['tick_s'], ha=tick_ha, va='top', fontsize=tick_fontsize, rotation=tick_rotation, color=row['color'])
            ax.text(row['acum'], row['bottom']-0.1, row['tick_e'], ha=tick_ha, va='top', fontsize=tick_fontsize, rotation=tick_rotation, color=row['color'])
        else:
            ax.text(row['acum'] - row['len']/2, row['bottom']-0.2, '{0}:{1}-{2}'.format(row['chrom'], row['tick_s'], row['tick_e']), ha='center', va='top', fontsize=tick_fontsize, color=row['color'])
            

    ax.set_xlim([0, sum_len])
    ax.set_ylim([-0.7, intervals-0.3])

    for i in ['bottom', 'top','right', 'left']:
        ax.spines[i].set_color('none')
        ax.spines[i].set_linewidth(0)

    ax.tick_params(bottom =False,top=False,left=False,right=False)
    ax.set_xticklabels("")
    ax.set_yticklabels("")
    #ax.spines["bottom"].set_color('black')
    #ax.spines["bottom"].set_linewidth(1)










