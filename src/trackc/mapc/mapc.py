import numpy as np
import pandas as pd
from matplotlib.colors import Colormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #,CenteredNorm, SymLogNorm,PowerNorm, Normalize  
from matplotlib.cm import get_cmap
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import seaborn as sns

#import sys
#import cooler
#import collections.abc as cabc
#import matplotlib as mpl
#from ..plotting._utils import panel_grid, ColorLike
#from matplotlib import rcParams, patheffects
#from copy import copy


fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)
fruitpunch2 = sns.blend_palette(['white', 'blue'], as_cmap=True)

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

def hex2rgb(value):
    # convert hex to rgb
    value = value.lstrip('#')
    lv = len(value)
    rgb = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
    return np.array([rgb[0]/255, rgb[1]/255, rgb[2]/255, 1])

def colorC(cname="RdBu_r", 
           bottom_color="#ffffff", 
           bad_color="white", 
           over_color="white", 
           under_color="white", 
           alpha=0):
    #cname = ["twilight_shifted", "jet", "RdBu_r", "RdGy_r", "BrBG_r", "hot_r", "Spectral_r"]
    #cname = ["terrain_r", "ocean", "gist_earth", "gist_stern_r", "tab20b", "twilight"]
    from matplotlib.colors import ListedColormap
    if isinstance(cname, str):
        cmap=plt.get_cmap(cname)
    else:
        cmap=cname

    cmap.set_bad(color=bad_color, alpha=alpha)
    cmap.set_over(color=over_color, alpha=alpha)
    cmap.set_under(color=under_color, alpha=alpha)

    if bottom_color==None:
        return cmap

    bottom_color = hex2rgb(bottom_color)
    newcolors = cmap(np.linspace(0, 1, 256))
    #white = np.array([1, 1, 1, 1])
    newcolors[0, :] = bottom_color
    newcmap = ListedColormap(newcolors)

    return newcmap

def plot_pcolormesh( 
    mat: np.ndarray,
    ax: Axes,
    norm=None,
    clim=None,
    cmap: Union[Colormap, str, None] = None,
    trans_XY: bool=False
    ):

    N = mat.shape[1]
    # Transformation matrix for rotating the heatmap.
    A = np.array([(y, x) for x in range(N, -1, -1) for y in range(N + 1)])
    t = np.array([[1,0.5], [-1,0.5]])
    A = np.dot(A, t)
  
    C = mat
    X = A[:, 1].reshape(N + 1, N + 1)
    Y = A[:, 0].reshape(N + 1, N + 1)

    if trans_XY==False:
        caxes = ax.pcolormesh(X, Y, np.flipud(C), cmap=cmap, edgecolor='none', norm=norm, clim=clim, snap=True, linewidth=.001)
    else:
        caxes = ax.pcolormesh(Y, X, np.flipud(C), cmap=cmap, edgecolor='none', norm=norm, clim=clim, snap=True, linewidth=.001)
    
    return caxes

def mapC_triview(
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
    trans_XY: bool = False,
    mapType: Union[str, None] = 'triangle',
    map_order: int = 0,
    symmetric: bool = False,
    k = 1
    ): 
    """\
    Draw triangle view of the C data

    Parameters
    mat:

    cmap
    mapType
        contact heatmap type, can be one of ``['square', 'triangle', 'rectangle']``,
        if select ``rectangle`` the ``mapHeight`` parameter should be set.
    mapHeight: int or None
        Parameter for ``mapType: rectangle'``, ``mapType: triangle'`` also can be set.
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
    if mapType=='square':
        im = ax.matshow(mat, cmap=cmap, norm=norm, clim=clim)
    else:
        if mapType=='triangle':
            im = plot_pcolormesh(mat=mat, ax=ax, norm=norm, clim=clim, cmap=cmap, trans_XY=trans_XY)
 
        if mapType=='rectangle':
            im = plot_pcolormesh(mat=mat, ax=ax, norm=norm, clim=clim, cmap=cmap, trans_XY=trans_XY)

    mapC_colorbar(ax, im, trans_XY, map_order, height, mapType)

    return im

def mapC_label(ax, label, trans_XY, map_order, mapType, height, fontsize= 10, color='k', Pair_mat=False):
    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim() 

    x, y = xmax, ymax
    verticalalignment = 'top'
    horizontalalignment = 'right'

    if mapType == "square":
        if map_order==0:
            x, y = xmax, ymax
        else:
            x, y = xmin, ymin
            verticalalignment = 'bottom'
            horizontalalignment ='left'
    else: # triangle, rectangle
        if Pair_mat==False:
            if trans_XY == True:
                ax.set_xlabel(label, fontsize=fontsize, color=color)
            else:
                ax.set_ylabel(label, fontsize=fontsize, color=color)
            
            return
        
        if trans_XY == True:
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

def mapC_colorbar(ax, im, trans_XY, map_order, height, mapType):
    x0, y0, x0_width, y0_height = 0, 1.015, 1, 0.015
    if trans_XY == False:
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
    if trans_XY:
        orientation = 'horizontal'
    
    cax = ax.inset_axes([x0, y0, x0_width, y0_height])
    cbar = plt.colorbar(im, ax=ax, cax=cax, orientation=orientation)
    if trans_XY:
        cbar.ax.xaxis.tick_top()
        cbar.ax.tick_params(axis='x', labelrotation=90)
    else:
        cbar.ax.yaxis.tick_right()
    #cbar.set_label(label, fontsize=8)

def paraPair(para):
    if isinstance(para, list)==False:
        return([para, para])
    else:
        if len(para)==1:
            para[1] = para[0]
    return para

def mapC(
    mat: Union[np.ndarray, None]  = None,
    mat2: Union[np.ndarray, None] = None,

    cmap: Union[Sequence[Colormap], Sequence[str], Colormap, str, None] = fruitpunch,
    label: Union[Sequence[str], str, None] = None,
    label_fontsize: Union[Sequence[int], int] = 10,
    label_color: Union[Sequence[str], str, None] = 'k',

    logdata: Union[Sequence[bool], bool] = False, 
    maxrange: Union[Sequence[float], float]=None,
    minrange: Union[Sequence[float], float]=None,
    trim_range: Union[Sequence[float], float]=0.99,

    ax: Optional[Axes] = None,
    #na_color: ColorLike = None,
    mapType: Union[str, None] = 'triangle',
    height: int = 0,
    trans_XY: bool = False,
    symmetric: bool = False,
    ax_on=True
    ):
    """\
    Draw triangle view of the C data

    Parameters
    mat:

    cmap
    mapType
        contact heatmap type, can be one of ``['square', 'triangle', 'rectangle']``,
        if select ``rectangle`` the ``mapHeight`` parameter should be set.
    mapHeight: int or None
        Parameter for ``mapType: rectangle'``, ``mapType: triangle'`` also can be set.
        it means the heapmap hight bin number
    symmetric: bool
        If there is one of ``mat`` and ``mat2`` para is None, 
        value ``True`` means the  symmetric heatmap

    """

    #cmap = copy(get_cmap(cmap))
    #cmap.set_bad(na_color)
    #cmap2 = copy(get_cmap(cmap2))
    #cmap2.set_bad(na_color)
    cmap = paraPair(cmap)
    label = paraPair(label)
    label_fontsize = paraPair(label_fontsize)
    label_color = paraPair(label_color)
    logdata = paraPair(logdata)
    maxrange = paraPair(maxrange)
    minrange = paraPair(minrange)
    trim_range = paraPair(trim_range)

    k = 1
    k2 = -1
    map_order = None
    map_order2 = None

    if isinstance(mat, np.ndarray)==False:
        k2 = 0
    if isinstance(mat2, np.ndarray)==False:
        k = 0

    if isinstance(mat, np.ndarray) and isinstance(mat2, np.ndarray):
        map_order = 0
        map_order2 = 1
 
    if isinstance(mat, np.ndarray):
        im = mapC_triview(
            mat=mat,
            ax=ax,
            cmap=cmap[0],
            label=label[0],
            logdata=logdata[0], 
            maxrange=maxrange[0],
            minrange = minrange[0],
            trim_range=trim_range[0],
            height=height,
            trans_XY=trans_XY,
            mapType=mapType,
            map_order=map_order,
            symmetric=symmetric,
            k = k
        )

    if isinstance(mat2, np.ndarray):
        im2 = mapC_triview(
            mat=mat2,
            ax=ax,
            cmap=cmap[1],
            label=label[1],
            logdata=logdata[1], 
            maxrange=maxrange[1],
            minrange = minrange[1],
            trim_range=trim_range[1],
            height=height,
            trans_XY=trans_XY,
            mapType=mapType,
            map_order=map_order2,
            symmetric=symmetric,
            k = k2
        )
   
    if mapType in ['triangle', 'rectangle']:
        if isinstance(mat2, np.ndarray)==False:
            if symmetric==False:
                if height > 0:
                    if trans_XY==False:
                        ax.set_ylim(0,height)
                        if mapType == 'rectangle':
                            ax.set_xlim(height, mat.shape[0]-height)
                    else:
                        ax.set_xlim(0,height)
                        if mapType == 'rectangle':
                            ax.set_ylim(height, mat.shape[0]-height)
                else:
                    if trans_XY==False:
                        ax.set_ylim(bottom=0)
                    else:
                        ax.set_xlim(left=0)
            else:
                if height > 0:
                    if trans_XY==False:
                        ax.set_ylim(-height, height)
                        if mapType == 'rectangle':
                            ax.set_xlim(height, mat.shape[0]-height)
                    else:
                        ax.set_xlim(-height, height)
                        if mapType == 'rectangle':
                            ax.set_ylim(height, mat.shape[0]-height)

        if isinstance(mat, np.ndarray)==False:
            if symmetric==False:
                if height > 0:
                    if trans_XY==False:
                        ax.set_ylim(-height, 0)
                        if mapType == 'rectangle':
                            ax.set_xlim(height, mat2.shape[0]-height)
                    else:
                        ax.set_xlim(-height,0)
                        if mapType == 'rectangle':
                            ax.set_ylim(height, mat2.shape[0]-height)
                else:
                    if trans_XY==False:
                        ax.set_ylim(top=0)
                    else:
                        ax.set_xlim(right=0)

            else:
                if height > 0:
                    if trans_XY==False:
                        ax.set_ylim(-height, height)
                        if mapType == 'rectangle':
                            ax.set_xlim(height, mat2.shape[0]-height)
                    else:
                        ax.set_xlim(-height, height)
                        if mapType == 'rectangle':
                            ax.set_ylim(height, mat2.shape[0]-height)

    Pair_mat = False
    if isinstance(mat, np.ndarray) and isinstance(mat2, np.ndarray):
        if height > 0:
            if trans_XY==False:
                ax.set_ylim(-height, height)
            else:
                ax.set_xlim(-height, height)
            
        if mapType == 'rectangle':
            if trans_XY==False:
                ax.set_xlim(height, mat2.shape[0]-height)
            else:
                ax.set_ylim(height, mat2.shape[0]-height)

        Pair_mat = True

    if isinstance(mat, np.ndarray):
        mapC_label(ax, label[0], trans_XY, 0, mapType, height, fontsize=label_fontsize[0], color=label_color[0], Pair_mat=Pair_mat)
    if isinstance(mat2, np.ndarray):
        mapC_label(ax, label[1], trans_XY, 1, mapType, height, fontsize=label_fontsize[1], color=label_color[1], Pair_mat=Pair_mat)

    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([]) 
    if ax_on==False:
        ax.axis('off')

def plot_contact_regions_line(ax, 
                              rl_regions, 
                              map_order=0, 
                              mapType='triangle', 
                              trans_XY=False, 
                              linewidth=1, 
                              linestyle='--', 
                              linecolor='k'):
    if trans_XY == True:
        print('trans_XY should be False')

    xA, yA, xB, yB = 0,0,0,0
    xA2, yA2, xB2, yB2 = 0,0,0,0
    rl_regions = rl_regions.query('chrom == chrom')

    cbins = rl_regions['cbins'].to_list()
    full_len = sum(cbins)
    v_sum = 0
    for i, v in enumerate(cbins):
        v_sum = v_sum+v

        if i == len(cbins)-1:
            break

        #v2 = cbins[i+1]

        if mapType == 'triangle':
            if map_order==0:
                xA, yA, xB, yB = v_sum/2, v_sum, v_sum, 0
                xA2, yA2, xB2, yB2 = v_sum, 0, v_sum+(full_len-v_sum)/2, full_len-v_sum
            else:
                xA, yA, xB, yB = v_sum/2, -v_sum, v_sum, 0
                xA2, yA2, xB2, yB2 = v_sum, 0, v_sum+(full_len-v_sum)/2, -(full_len-v_sum)
        
        if mapType == 'square':
            if map_order==0:
                xA, yA, xB, yB = v_sum-0.5, -0.5, v_sum-0.5, v_sum-0.5
                xA2, yA2, xB2, yB2 = v_sum-0.5, v_sum-0.5, full_len-0.5, v_sum-0.5
            else:
                xA, yA, xB, yB = -0.5, v_sum-0.5, v_sum-0.5, v_sum-0.5
                xA2, yA2, xB2, yB2 = v_sum-0.5, v_sum-0.5, v_sum-0.5, full_len-0.5

        ax.plot([xA, xB], [yA, yB], linestyle, linewidth=linewidth, color=linecolor)
        ax.plot([xA2, xB2], [yA2, yB2], linestyle, linewidth=linewidth, color=linecolor)



def plot_heatmap_triangle_xticks(ax, regin1_binN, regin2_binN, chrom1, start1, end1, chrom2, start2, end2, showXticks):
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

   
def colorbar_triangle(axm,im,ymax):
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


"""
def plot_chrom_arrow(ax, region_len, rev_region, regions, colors=my23colors[3:22]):
    ax.axis("off")
    ax.set_xlim(0,sum(region_len))
    ax.set_ylim(-2,2)
    
    start = 0
    for i, v in enumerate(region_len):
        if rev_region[i] == True:
            arrow = mpatches.FancyArrowPatch((start+v+1, 0), (start-1, 0), mutation_scale=30, color=colors[i])
            ax.add_patch(arrow) 
        else:
            arrow = mpatches.FancyArrowPatch((start-1, 0), (start+v+1, 0), mutation_scale=30, color=colors[i])
            ax.add_patch(arrow) 
            
        #ax.add_patch(arrow)    
        chrom, start_str, end_str = split_region(regions[i])
        #ax.text(start+region_len[i]/2, -1.5, chrom, horizontalalignment='center', verticalalignment='top', fontsize=10, color=colors[i])
        ax.text(start+region_len[i]/2, 0, chrom, horizontalalignment='center', verticalalignment='center', fontsize=10, color="black")
        
        ax.text(start, -1.5, start_str, horizontalalignment='left', verticalalignment='top', fontsize=9, color=colors[i], rotation=90)
        ax.text(start+region_len[i], -1.5, end_str, horizontalalignment='right', verticalalignment='top', fontsize=9, color=colors[i], rotation=90)
        start = start + v

"""    
    
        

