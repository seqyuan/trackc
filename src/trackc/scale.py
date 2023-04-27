from matplotlib.axes import Axes
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
from .mapc._getRegionsCmat import GenomeRegion

def scale_ticks(ticks, scale='bp', tick_fl='%0.2f'):
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
                regions: Union[Sequence[str], str, None] = None,
                chrom_fontsize: Union[int, None] = 10,
                scale_adjust: Union[str, None] = 'kb',
                tick_fl: Union[str, None] ='%0.2f',
                tick_fontsize: Union[int, None] = 8,
                tick_rotation: Union[int, None] = 0,
                ):
    """
    Parameters
    ----------
    
    scale_adjust
        ``str``: should be one of [kb, Mb, bp]
            adjust the scale unit to make it pretty
    tick_fl
        ``str``: The position label retains the number of decimal places                       

    """
    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    line_GenomeRegions = line_GenomeRegions.reset_index()

    chrom = line_GenomeRegions.loc[0, "chrom"]
    start = line_GenomeRegions.loc[0, "start"]
    end = line_GenomeRegions.loc[0, "end"]

    ax.set_xlim([start, end])

    xticks = ax.get_xticks()
    xtick_labels = xticks
    if scale_adjust == 'Mb':
        xtick_labels = xtick_labels/1000000
        xtick_labels = ["{0}".format(tick_fl % i) for i in xtick_labels]
        ax.set_xticks(xticks, xtick_labels, fontsize=tick_fontsize, rotation=tick_rotation)
        ax.spines['bottom'].set_position(('data', 0))
        ax.text(end-(abs(end-start)*0.05), -1, 'Mb', fontsize=chrom_fontsize)
    elif scale_adjust == 'kb':
        xtick_labels = xtick_labels/1000
        xtick_labels = ["{0}".format(tick_fl % i) for i in xtick_labels]
        ax.set_xticks(xticks, xtick_labels, fontsize=tick_fontsize, rotation=tick_rotation)
        ax.text(end-(abs(end-start)*0.05), -1, 'kb', fontsize=chrom_fontsize)
        ax.spines['bottom'].set_position(('data', 0))
    else:
        pass

    spines = ['top', 'right', 'left']
    for i in spines:
        ax.spines[i].set_visible(False)

    ax.set_yticks([])
    ax.set_yticklabels('')  
    ax.set_xlabel(chrom, fontsize=chrom_fontsize)
    ax.set_ylim([-1, 1])


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
    
    track_colors2= ['#EC6F5D', '#73C8DB', '#4EB4A2', '#6676A0', '#F5AE9B', '#9EA7C6', '#A7DBCC', '#E63634', '#9B806E']
    track_colors = ['#332488', '#347834', '#4AAB9A', '#88CCEE', '#DDCC77', '#CD6777', '#AA449A', '#892355', '#2572B2', '#D55F30', '#E6A03D']
    if colors != None:
        track_colors  = colors

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

    line_GenomeRegions['tick_s'] = scale_ticks(line_GenomeRegions['start'], scale=scale_adjust, tick_fl=tick_fl)
    line_GenomeRegions['tick_e'] = scale_ticks(line_GenomeRegions['end'], scale=scale_adjust, tick_fl=tick_fl)
    if scale_adjust in ['Mb', 'kb']:
        line_GenomeRegions.loc[line_GenomeRegions.index[-1], "tick_e"] = line_GenomeRegions.loc[line_GenomeRegions.index[-1], "tick_e"] + "(" + scale_adjust + ")"

    for i, row in line_GenomeRegions.iterrows():
        arrow_s = row['acum'] - row['len']
        dx = row['len']
        if row["isReverse"] == True:
            arrow_s = row['acum']
            dx = -1 * row['len']
            
        ax.arrow(arrow_s, row['bottom'], dx, 0, 
            overhang=1, width=0.01,
            head_width=0.25,
            head_length=sum_len/170,
            length_includes_head=True,
            color=row['color'],
            linewidth=1,
            alpha = alpha
            )
        
        tick_ha = 'center'
        if tick_rotation != 0:
            tick_ha = 'right'

        ax.text(row['acum'] - row['len']/2, row['bottom']+0.1, row['chrom'], ha='center', va='bottom', fontsize=tick_fontsize, color=row['color'])
        ax.text(row['acum'] - row['len'], row['bottom']-0.1, row['tick_s'], ha=tick_ha, va='top', fontsize=tick_fontsize, rotation=tick_rotation, color=row['color'])
        ax.text(row['acum'], row['bottom']-0.1, row['tick_e'], ha=tick_ha, va='top', fontsize=tick_fontsize, rotation=tick_rotation, color=row['color'])


    ax.set_xlim([0, sum_len])
    ax.set_ylim([-0.5, intervals-0.5])

    for i in ['bottom', 'top','right', 'left']:
        ax.spines[i].set_color('none')
        ax.spines[i].set_linewidth(0)

    ax.tick_params(bottom =False,top=False,left=False,right=False)
    ax.set_xticklabels("")
    ax.set_yticklabels("")
    #ax.spines["bottom"].set_color('black')
    #ax.spines["bottom"].set_linewidth(1)










