from matplotlib.axes import Axes
#import pyBigWig
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import pandas as pd
from trackc.tl._getRegionsCmat import GenomeRegion

def _make_multi_region_ax(ax, lineGenomeRegions):
    lineGenomeRegions['len'] = lineGenomeRegions['fetch_end']-lineGenomeRegions['fetch_start']
    lineGenomeRegions['ax_ratio'] = lineGenomeRegions['len']/lineGenomeRegions['len'].sum()
    lineGenomeRegions['ax_x'] = lineGenomeRegions['ax_ratio'].cumsum(axis=0) - lineGenomeRegions['ax_ratio']
    axs = [ax.inset_axes([row['ax_x'], 0, row['ax_ratio'], 1]) for i, row in lineGenomeRegions.iterrows()]
    for axi in axs:
        axi.axis('off')
    
    return axs

def bw_track(bw, 
             ax: Optional[Axes] = None,
             ylabel: Optional[str] = None,
             regions: Union[Sequence[str], str, None] = None, 
             binsize: Optional[int] = 50000,
             averagetype: Union[str, None] = 'mean',
             ymin: Optional[float] = None,
             ymax: Optional[float] = None,
             color: Union[Sequence[str], None] = '#827DBB',
             invert_y: Optional[bool] = False,
             label_rotation=0,
             label_fontsize: Optional[int] = 12,
             tick_fontsize: Optional[int] = 8,
             tick_fl: Optional[str] ='%0.2f', 
            ):
    """\
    Plot multi-regions bigwig signal tracks.
    
    Parameters
    ----------
    ax 
        ``cooler.Cooler``: cool format Hi-C matrix (https://github.com/open2c/cooler)
    ylabel
        ``str``: The ``'balance'`` parameters of ``coolMat.matrix(balance=False).fetch('chr6:119940450-123940450')``
    regions: bool, optional
        Force balancing weights to be interpreted as divisive (True) or
        multiplicativ

    binsize
        ``chrom region`` list: or ``chrom region`` or None. 
        The subset matrix row genome regions
        eg. ``"chr6:1000000-2000000"``, eg. ``["chr6:1000000-2000000", "chr3:5000000-4000000", "chr5"]``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means you want to get the reverse region contact matrix

    averagetype
        ``chrom region`` list: or ``chrom region`` or None. 
        The subset matrix col genome regions, default is ``None``, which means the sample region as ``row_regions``
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
    
    min_y = 0
    max_y = 0
    
    for i, row in line_GenomeRegions.iterrows():    
        bins = int(row['len']/binsize)
        plot_list = bw.stats(row['chrom'], row['fetch_start'], row['fetch_end'], type=averagetype, nBins=bins)
        plot_list = [0 if v is None else v  for v in plot_list]

        axs[i].bar(x=range(0, bins), height=plot_list, width=1, bottom=[0]*(bins),color=color[i],align="edge",edgecolor=color[i])    
        
        right, left = bins, 0
        if row['isReverse'] == True:
            left, right = bins, 0
        axs[i].set_xlim(left, right)
        
        if min_y < min(plot_list):
            min_y = min(plot_list)
            
        if max_y < max(plot_list):
            max_y = max(plot_list)
        
    if ymin == None:
        ymin = min_y
    if ymax == None:
        ymax = max_y
        
    if invert_y == True:
        ymin = max_y
        ymax = min_y

    for axi in axs:
        axi.set_ylim(ymin, ymax)
        
    ax.set_ylim(ymin, ymax)
    
    va = 'top'
    if invert_y == True:
        va='bottom'
    ax.text(0, ymax, " [{0}, {1}]".format(tick_fl % min_y, tick_fl % ymax), verticalalignment=va, fontsize=tick_fontsize)
    
    ax.set_ylabel(ylabel, fontsize=label_fontsize, rotation=label_rotation, horizontalalignment='right', verticalalignment='center')
     
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
    
def bw_compartment(compartment_bw, ax, chrom, start, end, ylabel, xticklabel=False, Acolor="#3271B2", Bcolor="#FBD23C", binsize=100000):
    chrsize = compartment_bw.chroms()[chrom]
    xbins = int(chrsize/binsize)
    if chrsize % binsize > 0:
        xbins = xbins + 1
        
    plot_list = compartment_bw.stats(chrom, 0, chrsize, nBins=xbins)
    mat=pd.DataFrame({"pc1":plot_list})
    mat["start"] = mat.index * binsize
    mat["width"] = binsize
    mat.loc[mat.index[-1],"length"] = chrsize % binsize
    mat.loc[mat[np.isnan(mat.pc1)].index, "pc1"] = 0
    
    plus = mat[mat['pc1']>0]
    minux = mat[mat['pc1']<=0]
 
    ax.bar(x=list(plus["start"]), height=plus['pc1'], width=plus["width"], bottom=[0]*(plus.shape[0]),color=Acolor,align="edge",edgecolor=Acolor,label="A")
    ax.bar(x=list(minux["start"]), height=minux['pc1'], width=minux["width"], bottom=[0]*(minux.shape[0]),color=Bcolor,align="edge",edgecolor=Bcolor,label="B")
    #ax.bar(0,height=0,color="#E27678",align="edge",edgecolor="#E27678",label="A2B")
    #ax.bar(0,height=0,color="#85AFBD",align="edge",edgecolor="#85AFBD",label="B2A")

    #ax.set_xlim(0, chrsize)
    ax.grid(False)
    ax.tick_params(bottom =False,top=False,left=True,right=False) #去掉tick线
    ax.spines['left'].set_color('k')
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.plot([0, chrsize], [0,0], '-', label='', linewidth=1, color='black', solid_capstyle='butt')
    #ax.set_yticklabels('')
    ax.set_ylabel(ylabel, fontsize=10, rotation='horizontal', horizontalalignment='right',verticalalignment='center')
    #ax.set_ylim([-1,1])
    #ax.set_yticks([-1, 1])
    ax.set_yticklabels([-1, 1], fontsize=8)
    if xticklabel==False:
        ax.set_xticklabels('')
    
    ax.set_xlim(start, end)

    
