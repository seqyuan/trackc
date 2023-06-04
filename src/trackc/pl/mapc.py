import numpy as np
import pandas as pd
from matplotlib.colors import Colormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #,CenteredNorm, SymLogNorm,PowerNorm, Normalize  
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
from trackc.pa import fruitpunch, fruitpunch2
import os
basedir = os.path.abspath(os.path.dirname(__file__))

ColorLike = Union[str, Tuple[float, ...]]

def _check_na_color(
    na_color: Optional[ColorLike], *, img: Optional[np.ndarray] = None
) -> ColorLike:
    if na_color is None:
        if img is not None:
            na_color = (0.0, 0.0, 0.0, 0.0)
        else:
            na_color = "lightgray"
    return na_color

def getData2Map(mat, maxrange=None, minrange=None, trim_range=0.99, inplace=False):
    mat = mat.astype(float)
    #if logdata: mat = np.log2(mat)   
    mat[mat == np.inf] = 0.0
    mat[mat == -np.inf] = 0.0
    mat[mat == 0.0] = np.nan    
    
    df = pd.DataFrame(np.ravel(mat))
    df = df[df != np.nan]
    
    if np.nanmax(mat)>0:
        if trim_range <1 and maxrange==None and minrange==None:
            xmaxrange = np.nanpercentile(abs(df), trim_range*100)
            xminrange = np.nanpercentile(abs(df), 100 - trim_range*100)
            print("no max min range")
        else:
            if maxrange==None:
                xmaxrange = abs(df).max()[0]
            else:
                xmaxrange = maxrange
            if minrange==None:
                xminrange = abs(df).min()[0]
            else:
                xminrange = minrange

        if inplace==True:
            mat[(mat<=xminrange) & (mat>0)] = xminrange    
            mat[(mat>=xmaxrange) & (mat>0)] = xmaxrange
            mat[(mat>= -xminrange) & (mat<0)] = -xminrange
            mat[(mat<= -xmaxrange) & (mat<0)] = -xmaxrange

        return mat, xmaxrange, xminrange
    else:
        print ("Warning: max data <0, no data to plot")


def _plot_pcolormesh( 
    mat: np.ndarray,
    ax: Axes,
    norm=None,
    clim=None,
    cmap: Union[Colormap, str, None] = None,
    trans_ax: bool=False
    ):

    N = mat.shape[1]
    # Transformation matrix for rotating the heatmap.
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[1,0.5], [-1,0.5]])
    A = np.dot(A, t)
  
    C = mat
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)

    if trans_ax==False:
        caxes = ax.pcolormesh(X, Y, np.flipud(C), cmap=cmap, edgecolor='none', norm=norm, clim=clim, snap=True, linewidth=.001)
    else:
        caxes = ax.pcolormesh(Y, X, np.flipud(C), cmap=cmap, edgecolor='none', norm=norm, clim=clim, snap=True, linewidth=.001)
    
    return caxes

def _mapC_triview(
    mat: np.ndarray,
    ax: Optional[Axes] = None,
    cmap: Union[Colormap, str, None] = fruitpunch,
    label: Union[str, None] = None,
    label_fontsize: int = 10,
    label_color: Union[str, None] = 'k',
    logdata: bool = False, 
    trim_range: float = 0.98,
    maxrange: float = None,
    minrange: float = None,
    height: int = 0,
    trans_ax: bool = False,
    map_type: Union[str, None] = 'triangle',
    map_order: int = 0,
    symmetric: bool = False,
    k = 1
    ): 
    """\
    Draw triangle view of the C data

    Parameters
    mat:

    cmap
    map_type
        contact heatmap type, can be one of ``['square', 'triangle', 'rectangle']``,
        if select ``rectangle`` the ``mapHeight`` parameter should be set.
    mapHeight: int or None
        Parameter for ``map_type: rectangle'``, ``map_type: triangle'`` also can be set.
        it means the heapmap hight bin number

    map_order: int
        0: np.triu()
        1: np.tril()
    k: int
        ``np.tril`` or ``np.triu`` ``k`` parameter

    """
    mat, maxrange, minrange = getData2Map(mat, maxrange=maxrange, minrange=minrange, trim_range=trim_range, inplace=False)
    print("maxrange:", maxrange ,"minrange:",minrange)

    if symmetric==False:
        if map_order==0:
            mat = np.triu(mat, k = k)
            mat[mat==0.0] = np.nan
        elif map_order==1:
            mat = np.tril(mat, k = k)
            mat[mat==0.0] = np.nan
        else:
            pass
    
    norm = None
    if logdata:
        norm=LogNorm()
    
    clim = (minrange, maxrange)
    if map_type=='square':
        im = ax.matshow(mat, cmap=cmap, norm=norm, clim=clim)
    else:
        if map_type=='triangle':
            im = _plot_pcolormesh(mat=mat, ax=ax, norm=norm, clim=clim, cmap=cmap, trans_ax=trans_ax)
 
        if map_type=='rectangle':
            im = _plot_pcolormesh(mat=mat, ax=ax, norm=norm, clim=clim, cmap=cmap, trans_ax=trans_ax)

    _mapC_colorbar(ax, im, trans_ax, map_order, height, map_type)

    return im

def _mapC_label(ax, label, trans_ax, map_order, map_type, height, fontsize= 10, color='k', Pair_mat=False):
    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim() 

    x, y = xmax, ymax
    verticalalignment = 'top'
    horizontalalignment = 'right'

    if map_type == "square":
        if map_order==0:
            x, y = xmax, ymax
        else:
            x, y = xmin, ymin
            verticalalignment = 'bottom'
            horizontalalignment ='left'
    else: # triangle, rectangle
        if Pair_mat==False:
            if trans_ax == True:
                ax.set_xlabel(label, fontsize=fontsize, color=color)
            else:
                ax.set_ylabel(label, fontsize=fontsize, color=color)
            
            return
        
        if trans_ax == True:
            y = ymax
            if map_order==0:
                x = xmax
            else:
                x = xmin
                horizontalalignment ='left'
        else:
            x = xmin
            horizontalalignment ='left'
            if map_order==1:
                y = ymin
                verticalalignment = 'bottom'

    bbox_props = dict(boxstyle="round", fc="w", alpha=0.9, pad=0, ec='none')

    ax.text(x, y, 
            label, 
            verticalalignment = verticalalignment,
            horizontalalignment = horizontalalignment,
            #transform = ax.transAxes,
            color =color, 
            fontsize = fontsize,
            bbox = bbox_props
            )

def _mapC_colorbar(ax, im, trans_ax, map_order, height, map_type):
    x0, y0, x0_width, y0_height = 0, 1.015, 1, 0.015
    if trans_ax == False:
        x0, y0, x0_width, y0_height = 1.015, 0, 0.015, 1
        if map_order in [0, 1]:
            y0_height = 0.45
            if height > 0 :
                x0_width = 0.012
            else:
                x0_width = 0.04
            if map_order==0:
                y0 = 0.55
            else:
                y0 = 0
        else:
            x0_width = 0.04
            if height > 0 :
                x0_width = 0.012
    else:
        if map_order in [0, 1]:
            x0_width = 0.45
            if height > 0 :
                y0_height = 0.012
            else:
                y0_height = 0.04
            if map_order==0:
                x0 = 0.55
            else:
                x0 = 0
        else:
            y0_height = 0.04
            if height > 0 :
                y0_height = 0.012

    orientation = 'vertical'
    if trans_ax:
        orientation = 'horizontal'
    
    cax = ax.inset_axes([x0, y0, x0_width, y0_height])
    cbar = plt.colorbar(im, ax=ax, cax=cax, orientation=orientation)
    if trans_ax:
        cbar.ax.xaxis.tick_top()
        cbar.ax.tick_params(axis='x', labelrotation=90)
    else:
        cbar.ax.yaxis.tick_right()
    #cbar.set_label(label, fontsize=8)

def _paraPair(para):
    if isinstance(para, list)==False:
        return([para, para])
    else:
        if len(para)==1:
            para[1] = para[0]
        else:
            return para
    return para

def mapC(
        ax: Optional[Axes] = None,
        mat: Union[np.ndarray, None]  = None,
        mat2: Union[np.ndarray, None] = None,
        
        cmap: Union[Sequence[Colormap], Sequence[str], Colormap, str, None] = [fruitpunch, 'YlOrRd'],
        label: Union[Sequence[str], str, None] = None,
        label_fontsize: Union[Sequence[int], int] = 10,
        label_color: Union[Sequence[str], str, None] = 'k',
        
        logdata: Union[Sequence[bool], bool] = False, 
        maxrange: Union[Sequence[float], float]=None,
        minrange: Union[Sequence[float], float]=None,
        trim_range: Union[Sequence[float], float]=0.98,
        
        map_type: Union[str, None] = 'triangle',
        height: int = 0,
        trans_ax: bool = False,
        symmetric: bool = False,
        ax_on: bool =True,
        aspect: Union[str, float]='auto'
        ):
    """\
    Plot contact map, support for multiple or reverse genome regions.
    This function implements the plot method for `np.ndarray`,
    which could get from trackc.tl.extractCisContact or trackc.tl.extractContactRegions
    By default, the trim_range value is fixed so that the 98th percentile (resp. 2th percentile) of each
    interaction matrix is discarded. It therefore allow to remove the extreme values from the matrix,  
    mat or mat2 is plotted independently
    If the maxrange parameter is set, data higher that this threshold will be fixed to the maxrange value.
    
    cmap, label, label_fontsize, label_color, logdata, minrange, maxrange, trim_range, those parameters can be set as a list, 
    mat and mat2 will set to the first two values. If those parameters are one single value or the length of list is one,
    then mat and mat2 both set the same value
   
    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
    mat: `np.ndarray`
        matrix for plot upper or right of heatmap
    mat2: `np.ndarray`
        matrix for plot bottom or left of heatmap
    cmap: `str` | `matplotlib.colors.Colormap` | `list`
        colormap for continuous annotations, if set as list, mat and mat2 will set to the first two values
    label: `str`
        the title of the track, will show on the left
    label_fontsize: `int`
        the label text fontsize
    label_color: `int`
        the label text color
    logdata: `bool` | `bool list`
        do you want to log the data before plotting the heatmap
    minrange: `float` | `float list`
        the minimum range of values used to define the color palette
    maxrange: `float` | `float list` 
        the maximum range of values used to define the color palette
    trim_range: `float` | `float list` 
        remove the extreme values by trimming the counts.[0,1]
        define the maxrange and minrange values using the percentile of the interaction matrix
    map_type: `str`
        optional is ['square', 'triangle', 'rectangle'], default is square
        triangle and rectangle default is flip the image 45 degrees to the left.
        for rectangle type, the corresponding length of ``height`` will be truncated 
        from both ends of the input matrix.
    height: `int`
        if map_type is one of ['triangle', 'rectangle'], `height` indicates the longest interaction bin interval you want to show
    trans_ax: `bool`
        whether flip the image 45 degrees to the right
    symmetric: `bool`
        whether to display a symmetrical heatmap when only one of mat and mat2 is set
    ax_on: `bool`
        whether show the spines
    aspect: `str` | `float`
        optional is 'auto' or 1
        length-width ratio of heatmap
    
    Example
    -------
    >>> import trackc as tc
    >>> import cooler
    >>> BxPC3 = cooler.Cooler('./BxPC3.chr18.mcool::/resolutions/25000')
    >>> neo_domain_regions = ['18:47950000-48280000', '18:75280000-74850000']
    >>> tumor_zoom = tc.tl.extractContactRegions(clr=BxPC3, row_regions=neo_domain_regions)
    >>> ten = tc.tenon(width=6, height=1)
    >>> ten.add(pos='bottom', height=0.7, hspace=0.05)
    >>> tc.pl.mapC(ax=ten.axs(0), mat=tumor_zoom.cmat, map_type='triangle',
            maxrange=200, minrange=10, label='tumor res=25k', ax_on=False, height=40)
    >>> tc.savefig('trackc_mapc.pdf')
    """

    cmap = _paraPair(cmap)
    label = _paraPair(label)
    label_fontsize = _paraPair(label_fontsize)
    label_color = _paraPair(label_color)
    logdata = _paraPair(logdata)
    maxrange = _paraPair(maxrange)
    minrange = _paraPair(minrange)
    trim_range = _paraPair(trim_range)

    k = 1
    k2 = -1

    if isinstance(mat, np.ndarray)==False:
        k2 = 0
    if isinstance(mat2, np.ndarray)==False:
        k = 0

    #ax2 = ax
    Pair_mat = False
    ax2 = ax.inset_axes([0, 0, 1, 1], facecolor='none')
    ax2.set_zorder(1)
    ax2.set_xticklabels([])
    ax2.set_xticks([])
    ax2.set_yticklabels([])
    ax2.set_yticks([]) 


    if isinstance(mat, np.ndarray) and isinstance(mat2, np.ndarray):
        #ax2 = ax.inset_axes([0, 0, 1, 1], facecolor='none')
        #ax2.set_zorder(-1)
        symmetric = False
        Pair_mat = True

    if isinstance(mat, np.ndarray):
        im = _mapC_triview(
            mat=mat,
            ax=ax,
            cmap=cmap[0],
            label=label[0],
            logdata=logdata[0], 
            maxrange=maxrange[0],
            minrange = minrange[0],
            trim_range=trim_range[0],
            height=height,
            trans_ax=trans_ax,
            map_type=map_type,
            map_order=0,
            symmetric=symmetric,
            k = k
        )

    if isinstance(mat2, np.ndarray):
        #ax2 = ax.inset_axes([0, 0, 1, 1], facecolor='none')
        #ax2.set_zorder(-1)
        
        im2 = _mapC_triview(
            mat=mat2,
            ax=ax2,
            cmap=cmap[1],
            label=label[1],
            logdata=logdata[1], 
            maxrange=maxrange[1],
            minrange = minrange[1],
            trim_range=trim_range[1],
            height=height,
            trans_ax=trans_ax,
            map_type=map_type,
            map_order=1,
            symmetric=symmetric,
            k = k2
        )
        
    xylim_conf = pd.read_table(os.path.join(basedir,'mapc_xylim.txt'), header=0)

    if isinstance(mat, np.ndarray) and isinstance(mat2, np.ndarray):
        set_xylim(ax, xylim_conf, map_type, 0, symmetric, trans_ax, Pair_mat, mat, height)
        set_xylim(ax2, xylim_conf, map_type, 1, symmetric, trans_ax, Pair_mat, mat2, height)

        ax.set_aspect(aspect)
        ax2.set_aspect(aspect)
        
    elif isinstance(mat, np.ndarray):
        set_xylim(ax, xylim_conf, map_type, 0, symmetric, trans_ax, Pair_mat, mat, height)
        ax.set_aspect(aspect)
    elif isinstance(mat2, np.ndarray):
        set_xylim(ax, xylim_conf, map_type, 1, symmetric, trans_ax, Pair_mat, mat2, height)
        set_xylim(ax2, xylim_conf, map_type, 1, symmetric, trans_ax, Pair_mat, mat2, height)
        ax.set_aspect(aspect)
        ax2.set_aspect(aspect)
    else:
        pass

    if ax_on==False:
        ax.axis('off')
        ax2.axis('off')

    if isinstance(mat, np.ndarray):
        _mapC_label(ax, label[0], trans_ax, 0, map_type, height, fontsize=label_fontsize[0], color=label_color[0], Pair_mat=Pair_mat)
    if isinstance(mat2, np.ndarray):
        _mapC_label(ax2, label[1], trans_ax, 1, map_type, height, fontsize=label_fontsize[1], color=label_color[1], Pair_mat=Pair_mat)

def set_xylim(axs, conf, map_type, map_order, symmetric, trans_ax, pair_mat, matx, height):
    conf = conf[(conf['map_type']==map_type) & (conf['map_order']==map_order) & (conf['symmetric']==symmetric) & (conf['trans_ax']==trans_ax) & (conf['pair_mat']==pair_mat)]
    if height == 0:
        height = matx.shape[0]
    
    conf_values = {'mat_w': matx.shape[1], 'mat_h': matx.shape[0], 'height': height, 'height_m': -height, 
                   'mat_w_height': matx.shape[1]-height, '0':0, '-0.5':-0.5, 'mat_w_0.5m':matx.shape[1]-0.5, 'mat_h_0.5m':matx.shape[0]-0.5}
    
    if conf.shape[0] == 0:
        print('type not set xylim')
    else:
        ix = conf.index[0]
        xmin = conf_values[conf.loc[ix, 'xlim_left']]
        xmax = conf_values[conf.loc[ix, 'xlim_right']]
        ymin = conf_values[conf.loc[ix, 'ylim_bottom']]
        ymax = conf_values[conf.loc[ix, 'ylim_top']]
        axs.set_xlim([xmin, xmax])
        axs.set_ylim([ymin, ymax])
    
    axs.set_xticklabels([])
    axs.set_xticks([])
    axs.set_yticklabels([])
    axs.set_yticks([]) 
    
    

def _plot_heatmap_triangle_xticks(ax, regin1_binN, regin2_binN, chrom1, start1, end1, chrom2, start2, end2, showXticks):
    #ax.set_xticks([])
    ax.set_xticks([0, regin1_binN, regin1_binN + regin2_binN])
    ax.set_xticklabels([start1, str(end1) + " " + str(start2), end2])
    
    ax.set_yticklabels([])
    ax.set_yticks([])
    
    ax.xaxis.tick_top()
    ax.spines['top'].set_color('k')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')

   
def _colorbar_triangle(axm,im,ymax):
    height="2%"
    width="10%"
    if ymax!=None:
        height="6%"
    axins1 = inset_axes(axm, width=width, height=height, loc=4, bbox_to_anchor=(-0.1, 0.7, 1.2, 2.9), bbox_transform=axm.transAxes,)
    cbar=plt.colorbar(im, cax=axins1, orientation='horizontal')
    cbar.ax.xaxis.tick_top()
    cbar.ax.spines['top'].set_color('none')
    cbar.ax.spines['right'].set_color('none')
    cbar.ax.spines['bottom'].set_color('none')
    cbar.ax.spines['left'].set_color('none')



        

