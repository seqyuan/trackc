from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
import pandas as pd
import numpy as np
from ..tools._getRegionsCmat import GenomeRegion
from .bigwig import make_multi_region_ax


def two_degree_bc(x_l=10, x_r=90, y_lr=0, y2=10, dots_num=100): 
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

def plot_loop(ax, loop_df, color, max_extend, invert_y, start, end, left_anchor, right_anchor):
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
        
        xt, yt = two_degree_bc(x_l=row[left_anchor], x_r=row[right_anchor], y_lr=0, y2=top, dots_num=100)
    
        ax.plot(xt, yt, color=color, linewidth=0.5, solid_capstyle='butt')
        if max(yt) > top_y:
            top_y = max(yt)
            
    ax.set_xlim(start, end)
    if invert_y==True:
        ax.set_ylim(0.5, 0)
    else:
        ax.set_ylim(0, 0.5)

def links_track(links_df: pd.DataFrame, 
              ax: Optional[Axes] = None,
              ylabel: Optional[str] = None,
              regions: Union[Sequence[str], str, None] = None, 
              links_type: Union[str, None] = 'loop',
              color: Union[Sequence[str], None] = '#66AC84',
              maxrange: Union[int, None] = 3000000,
              invert_y: Optional[bool] = False,
              anchor: Union[str, None] = 'inside',
              label_rotation: Union[int, None] = 0,
              label_fontsize: Optional[int] = 12,
              ):
    """\
    Plot multi-regions loop or TAD links.
    
    Parameters
    ----------
    loop_bed
        ``pd.DataFrame``: the file format for links is (tab separated), column names:
            chr1 x1 x2 chr2 y1 y2 (score ...)
        The score field is optional
        The fields after the y2 score will be ignored, score column is useless right now
        for example:
            chr1 100 200 chr1 250 300 1

    ax 
        ``cooler.Cooler``: cool format Hi-C matrix (https://github.com/open2c/cooler)
    ylabel
        ``str``: The ``'balance'`` parameters of ``coolMat.matrix(balance=False).fetch('chr6:119940450-123940450')``
    regions:
        bool, optional
        Force balancing weights to be interpreted as divisive (True) or
        multiplicativ
    links_type: 
        ``str``: links type, either 'loop' or 'triangle'
    
    
    anchor:
        ``str``: links type, either 'loop' or 'triangle'
            the link coordinates type for loop and triangle, by default: inside
                inside: link x2 and y1
                mid: link the middle of x1 and x2 and the middle of y1 and y2, (x1+x2)/2 and (y1+y2)/2
                outside: link x1 and y2
    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    axs = make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    if isinstance(color, list)==False:
        color = [color]
    if len(color) < line_GenomeRegions.shape[0]:
        repeat_times = (line_GenomeRegions.shape[0] + len(color) - 1) // len(color)
        color = (color * repeat_times)[:line_GenomeRegions.shape[0]]
    
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

    links_df = links_df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']]
    links_df['chr1'] = links_df['chr1'].astype(str)
    links_df['chr2'] = links_df['chr2'].astype(str)

    for ix, row in line_GenomeRegions.iterrows(): 
        loop_bed_plot = links_df[links_df['chr1']==row['chrom']]
        

        loop_bed_plot['length'] = loop_bed_plot["y2"] - loop_bed_plot["x1"]
        
        loop_bed_plot = loop_bed_plot[((loop_bed_plot['x1'] >= row['fetch_start']) & (loop_bed_plot['x1'] <= row['fetch_end'])) | \
                                ((loop_bed_plot['x2'] >= row['fetch_start']) & (loop_bed_plot['x2'] <= row['fetch_end']) |\
                                (loop_bed_plot['y1'] >= row['fetch_start']) & (loop_bed_plot['y1'] <= row['fetch_end'])) | \
                                    ((loop_bed_plot['y2'] >= row['fetch_start']) & (loop_bed_plot['y2'] <= row['fetch_end']))]

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
        if links_type == 'loop':
            plot_loop(axs[ix], links_df_list[ix], color[ix], max_extend, invert_y, row['start'], row['end'], left_anchor, right_anchor)
        #if links_type == 'loop':
            #plot_triangle()

    ax.set_ylabel(ylabel, fontsize=label_fontsize, rotation=label_rotation, horizontalalignment='right',verticalalignment='center')
    
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

