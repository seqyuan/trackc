from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import Colormap
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
from trackc.tl._getRegionsCmat import GenomeRegion
from .links import _plot_loop_arc
from .bigwig import _make_multi_region_ax

def bed_track(ax: Optional[Axes] = None,
              bed: Union[pd.DataFrame, str, None]  = None,
              regions: Union[Sequence[str], str, None] = None,
              track_style: Union[str, None] = 'bar',
              color: Union[Sequence[str], None] = 'tab:blue',
              cmap: Union[Colormap, str, None] = None,
              intervals: Union[int, None] = 1,
              #show_names: Union[bool, None] = False,
              alpha: Union[float, None] = 1,
              invert_y: Optional[bool] = False,
              label: Union[str, None] = None,
              label_fontsize: Union[int, None] = 12,
              label_rotation: Union[int, None] = 0,
              ymin: Optional[float] = None,
              ymax: Optional[float] = None,
              tick_fontsize: Optional[int] = 8,
              tick_fl: Optional[str] ='%0.2f',
              score_label_size: Union[int, None] = 7,
              ):
    """\
    Plot bed track, support for multiple or reverse genome regions.
    support bed3 and bed5, the fields after the column5 will be ignored, 
        should be sorted py chromStart if ``track_style`` is `line`

    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    bed: `pd.DataFrame` | `str`
        If ``bed`` if a filepath, the file should have no headers
        Here is bed formats:
        column1: chrom
            The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        column2: chromStart
            The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        column3: chromEnd
            The ending position of the feature in the chromosome or scaffold.
        column4: name
            Defines the name of the BED line. Either "." (=no name), or other string
        column5: score
            Defines the name of the BED line. if nessasary, can be set as ``.``,
            if track_type/style is one of bar/line
        column6: strand 
            Defines the strand. Either "." (=no strand) or "+" or "-".

        if the input data have 4 columns, then default input is bedGraph format: [chrom start end score]

    regions: `str` | `str list`
        The genome regions to plot
        e.g. ``"chr6:1000000-2000000"`` or ``["chr6:1000000-2000000", "chr3:5000000-4000000"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means you want to get the reverse region
    track_style: `str`
        bed blocks style,  opions in ['line', 'bar', 'link', 'tri', 'rec']
    color: `str` or `list`
        the color of line/triangle/rectangle, if color is color list, the block will set by regions
    cmap: `str` | `matplotlib.colors.Colormap`
        the colormap of the plot except track_style:line
    intervals: 


    intervals
        ``int``: if track_style is one of [tri, rec], the row number distribution for triangle or rectangle blocks
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

    if isinstance(cmap, list)==False:
        cmap = [cmap]
    if len(cmap) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(cmap) - 1) // len(cmap)
        cmap = (cmap * repeat_times)[:line_GenomeRegions.shape[0]]

    ax.set_ylabel(label, fontsize=label_fontsize, rotation=label_rotation, ha='right', va='center')
    spines = ['top', 'right', 'left']
    for i in spines:
        ax.spines[i].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels('')
    ax.set_yticks([])
    ax.set_yticklabels('')  

    if isinstance(bed, str)==True:
        bed=pd.read_table(bed, sep="\t", header=None)

    score_label = ''
    if bed.shape[1] == 3:
        bed.columns = ['chrom', 'start', 'end']
    if bed.shape[1] == 4:
        bed.columns = ['chrom', 'start', 'end', 'score']
        if bed['score'].dtypes=='int64':
            print('Input data is bedGraph format')
    if bed.shape[1] == 5:
        score_label = bed.columns[4]
        bed.columns = ['chrom', 'start', 'end', 'name', 'score']
    if bed.shape[1] == 6:
        score_label = bed.columns[4]
        bed.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed['chrom'] = bed['chrom'].astype(str)

    min_y = None
    max_y = None
    max_len = 0
    
    for ix, row in line_GenomeRegions.iterrows(): 
        bed2plot = bed[(bed['chrom']==row['chrom']) & (bed['end']>=row['fetch_start']) & (bed['start']<=row['fetch_end'])]
        if track_style in ['line', 'bar'] or bed.shape[1]>=5:
            if min_y==None:
                min_y = bed2plot['score'].min(skipna=True, numeric_only=True)
            else:
                if min_y < bed2plot['score'].min(skipna=True, numeric_only=True):
                    min_y = bed2plot['score'].min(skipna=True, numeric_only=True)
        
            if max_y==None:
                max_y = bed2plot['score'].max(skipna=True, numeric_only=True)
            else:
                if max_y > bed2plot['score'].max(skipna=True, numeric_only=True):
                    max_y = bed2plot['score'].max(skipna=True, numeric_only=True)

        maxlength = (bed2plot['end']-bed2plot['start']).max(skipna=True, numeric_only=True)
        if max_len < maxlength:
            max_len = maxlength

    if ymin == None:
        ymin = min_y
    if ymax == None:
        ymax = max_y

    
    for ix, row in line_GenomeRegions.iterrows(): 
        bed2plot = bed[(bed['chrom']==row['chrom']) & (bed['end']>=row['fetch_start']) & (bed['start']<=row['fetch_end'])].copy()
        if bed2plot.shape[0] == 0:
            continue
        
        if track_style == "line":
            _plot_bed_bar_l(axs[ix], bed2plot, row['fetch_start'], row['fetch_end'], needReverse=row['isReverse'], style='line', color=color[ix], alpha=alpha)    
        
        if track_style == "bar":
            _plot_bed_bar_l(axs[ix], bed2plot, row['fetch_start'], row['fetch_end'], needReverse=row['isReverse'], style='bar', color=color[ix], alpha=alpha)
        
        if track_style == "rec":
            _plot_bed_rec(ax, axs[ix], bed2plot, row['fetch_start'], row['fetch_end'], 
                         needReverse=row['isReverse'], color=color[ix], cname=cmap[ix], 
                         alpha=alpha, min=ymin, max=ymax, score_label=score_label, intervals=intervals, score_label_size=score_label_size)
        
        if track_style == "tri":
            _plot_bed_tri(ax, axs[ix], bed2plot, row['fetch_start'], row['fetch_end'], 
                         needReverse=row['isReverse'], color=color[ix], cname=cmap[ix], 
                         alpha=alpha, min=ymin, max=ymax, score_label=score_label, score_label_size=score_label_size)
            
            if invert_y:
                axs[ix].set_ylim(max_len/2, 0)
            else:
                axs[ix].set_ylim(0, max_len/2)
        
        if track_style == "link":
            _plot_bed_link(ax=ax, bed=bed2plot, 
                           start=row['fetch_start'],end=row['fetch_end'], 
                           needReverse=row['isReverse'], invert_y=invert_y, 
                           color=color[ix], cmap=cmap[ix], alpha=alpha)

    if track_style in ['line', 'bar']:
        for axi in axs:
            if invert_y:
                axi.set_ylim(ymax, ymin)
            else:
                axi.set_ylim(ymin, ymax)

        if invert_y:
            ax.set_ylim(ymax, ymin)
        else:
            ax.set_ylim(ymin, ymax)
            
        text_label_y_pos = ymax
        if invert_y:
            text_label_y_pos = ymin
            ax.text(0, ymax, " [{1}, {0}]".format(tick_fl % ymin, tick_fl % ymax), va='bottom', fontsize=tick_fontsize)
        else:
            ax.text(0, text_label_y_pos, " [{0}, {1}]".format(tick_fl % ymin, tick_fl % ymax), va='top', fontsize=tick_fontsize)


def _make_tri_data(start, end):
    data = np.array([[start, 0], [end, 0], [start+(end-start)/2, (end-start)/2]])
    return data

def _plot_bed_tri(mainAX, ax, bed, start, end, needReverse, color, cname, alpha, min, max, score_label=None, score_label_size=8):
    colors = color
    norm = None
    if cname != None:
        if isinstance(cname, str):
            map_vir=cm.get_cmap(cname)
        if isinstance(cname, Colormap):
            map_vir=cname
        
        # 因为 y 大到一定程度超过临界数值后颜色就会饱和不变(不使用循环colormap)。
        norm = plt.Normalize(min, max)
        # matplotlib.colors.Normalize 对象，可以作为参数传入到绘图方法里
        # 也可给其传入数值直接计算归一化的结果
        norm_y = norm(bed['score'])
        colors = map_vir(norm_y)
        
    patches = []
    bed = bed.reset_index()
    #print(colors[0])
    #print(colors[1])

    for i, row in bed.iterrows():
        polygon = Polygon(_make_tri_data(row['start'], row['end']), True, color=colors)
        patches.append(polygon)
    
    p = PatchCollection(patches, alpha=alpha, match_original=True)
    ax.add_collection(p)
    
    xlim_s = start
    xlim_e = end
    if needReverse == True:
         xlim_s = end
         xlim_e = start

    ax.set_xlim(xlim_s, xlim_e)
    if cname != None:
        sm = cm.ScalarMappable(norm=norm, cmap=map_vir)
        cax = mainAX.inset_axes([1.01, 0, 0.01, 1])
        cb = plt.colorbar(sm, ax=mainAX, cax=cax, label=score_label)
        cb.set_label(score_label, fontsize=score_label_size)

def _plot_bed_rec(mainAX, ax, bed, start, end, needReverse, color, cname, alpha, min, max, score_label=None, intervals=1, score_label_size=8):
    colors = color
    if cname != None:
        if isinstance(cname, str):
            map_vir=cm.get_cmap(cname)
        if isinstance(cname, Colormap):
            map_vir=cname
        
        # 因为 y 大到一定程度超过临界数值后颜色就会饱和不变(不使用循环colormap)。
        norm = plt.Normalize(min, max)
        # matplotlib.colors.Normalize 对象，可以作为参数传入到绘图方法里
        # 也可给其传入数值直接计算归一化的结果
        norm_y = norm(bed['score'])
        colors = map_vir(norm_y)
    bed.loc[:,'length'] = abs(bed['end'] - bed['start'])
    bottom = list(range(intervals))  
    # broadcast
    if intervals < bed.shape[0]:
        repeat_times = (bed.shape[0] + intervals - 1) // intervals
        bottom = (bottom * repeat_times)[:bed.shape[0]]
    else:
        bottom = bottom[:bed.shape[0]]
    bed.loc[:, 'bottom'] = bottom

    ax.bar(x='start', width='length', height=1, bottom='bottom', align='edge', color= colors, edgecolor=None, data=bed, alpha=alpha)
    xlim_s = start
    xlim_e = end
    if needReverse == True:
         xlim_s = end
         xlim_e = start
    ax.set_xlim(xlim_s, xlim_e)
    if cname != None:
        sm = cm.ScalarMappable(norm=norm, cmap=map_vir)
        cax = mainAX.inset_axes([1.01, 0, 0.01, 1])
        cb = plt.colorbar(sm, ax=mainAX, cax=cax, label=score_label)
        cb.set_label(score_label, fontsize=score_label_size)


def _plot_bed_bar_l(ax, bed, start, end, needReverse, style='bar', color='tab:blue', alpha=1):
    if style == 'bar':
        ax.bar(x=bed['start'], width=bed['end']-bed['start'], height=bed['score'], color=color, alpha=alpha, align='edge')
    
    if style == 'line':
        ax.plot(bed['start'], bed['score'], color=color, alpha=alpha, solid_capstyle='butt')
    
    xlim_s = start
    xlim_e = end
    if needReverse == True:
         xlim_s = end
         xlim_e = start
    ax.set_xlim(xlim_s, xlim_e)
    
    for i in ['top','right', 'left']:
        ax.spines[i].set_color('none')
        ax.spines[i].set_linewidth(0)
    ax.spines["bottom"].set_color('black')
    ax.spines["bottom"].set_linewidth(1)
    #ax.tick_params(bottom =True,top=False,left=False,right=False)
    #ax.set_xticklabels("")
    #ax.set_yticklabels("")

def _plot_bed_link(ax, bed, start, end, needReverse=False, invert_y=False, color='tab:blue', cmap='RdBu', alpha=1):
    bed = bed[bed.columns[[0,1,2]]]
    bed.columns = ['chrom', 'start', 'end']
    bed['length'] = bed['end'] - bed['start']
    max_extend = 0
    if bed.shape[0] >0:
        if max_extend < max(bed['length']):
            max_extend = max(bed['length'])

    _plot_loop_arc(ax, bed, color, max_extend, invert_y, start, end, 'start', 'end')
    
    if needReverse == True:
        ax.set_xlim(end, start)
    else:
        ax.set_xlim(start, end)
    
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