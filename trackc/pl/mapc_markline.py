from trackc.tl._getRegionsCmat import GenomeRegion
from .Virtual4C import _get_pets
from matplotlib.axes import Axes
from typing import Union, Sequence, Optional
import pandas as pd

def _mark_contact_regions_square(ax, 
                                 markRegions, 
                                 plot_region_df, 
                                 map_order, 
                                 symmetric, 
                                 linewidth, 
                                 linestyle, 
                                 linecolor, 
                                 only_cis
                                 ):
    for i, reg in enumerate(markRegions):
        reg_pos = plot_region_df[plot_region_df[reg]==1]
        x0 = reg_pos.index[0] - 0.5
        x1 = reg_pos.index[-1] + 0.5
        
        """
        ---> ln1
        |ln3 |ln2
        ---->ln4 
        
        (pos1)   _____ (pos4)
                |     |
                |     |
        (pos2)   _____ (pos3)
        """
        ln1, = ax.plot((x0, x1), (x0, x0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        ln2, = ax.plot((x1, x1), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        ln3, = ax.plot((x0, x0), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        ln4, = ax.plot((x0, x1), (x1, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')

        if symmetric == False:
            if map_order == 0:
                ln3.remove()
                ln4.remove()
            else:
                ln1.remove()
                ln2.remove()

        if only_cis==False:
            for region in markRegions[i+1:]:
                trans_reg_pos = plot_region_df[plot_region_df[region]==1]
                y0 = trans_reg_pos.index[0]
                y1 = trans_reg_pos.index[-1]
                ln5, = ax.plot((x0, x1), (y0, y0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln6, = ax.plot((x0, x1), (y1, y1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                ln7, = ax.plot((y0, y0), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln8, = ax.plot((y1, y1), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                if symmetric == False:
                    if map_order == 1:
                        ln7.remove()
                        ln8.remove()
                    if map_order == 0:
                        ln5.remove()
                        ln6.remove()
            for region in markRegions[:i]:
                trans_reg_pos = plot_region_df[plot_region_df[region]==1]
                y0 = trans_reg_pos.index[0]
                y1 = trans_reg_pos.index[-1]
                ln5, = ax.plot((y0, y0), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln6, = ax.plot((y1, y1), (x0, x1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                ln7, = ax.plot((x0, x1), (y0, y0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln8, = ax.plot((x0, x1), (y1, y1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                if symmetric == False:
                    if map_order == 1:
                        ln7.remove()
                        ln8.remove()
                    if map_order == 0:
                        ln5.remove()
                        ln6.remove()
    

def _mark_contact_regions_triangle(ax, 
                                   markRegions, 
                                   plot_region_df, 
                                   map_order, 
                                   symmetric, 
                                   linewidth, 
                                   linestyle, 
                                   linecolor, 
                                   only_cis
                                   ):
    for i, reg in enumerate(markRegions):
        reg_pos = plot_region_df[plot_region_df[reg]==1]

        x0 = reg_pos.index[0]
        x1 = reg_pos.index[-1]
        ln1, = ax.plot((x0, x0+(x1-x0)/2), (0, x1-x0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        ln2, = ax.plot((x0+(x1-x0)/2, x1), (x1-x0, 0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        
        ln3, = ax.plot((x0, x0+(x1-x0)/2), (0, (x1-x0)*-1), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        ln4, = ax.plot((x0+(x1-x0)/2, x1), ((x1-x0)*-1, 0), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
        if symmetric == False:
            if map_order == 0:
                ln3.remove()
                ln4.remove()
            else:
                ln1.remove()
                ln2.remove()

        if only_cis==False:
            for region in markRegions[i+1:]:
                trans_reg_pos = plot_region_df[plot_region_df[region]==1]
                y0 = trans_reg_pos.index[0]
                y1 = trans_reg_pos.index[-1]
                pos_botom = (x1+(y0-x1)/2, y0-x1)
                pos_top = (x0+(y1-x0)/2, y1-x0)
                pos_l = (x0+(y0-x0)/2, y0-x0)
                pos_r = (x1+(y1-x1)/2, y1-x1)

                ln5, = ax.plot((pos_botom[0], pos_r[0]), (pos_botom[1], pos_r[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln6, = ax.plot((pos_top[0], pos_l[0]), (pos_top[1], pos_l[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln7, = ax.plot((pos_botom[0], pos_r[0]), (-pos_botom[1], -pos_r[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln8, = ax.plot((pos_top[0], pos_l[0]), (-pos_top[1], -pos_l[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                ln5_2, = ax.plot((pos_botom[0], pos_l[0]), (pos_botom[1], pos_l[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln6_2, = ax.plot((pos_top[0], pos_r[0]), (pos_top[1], pos_r[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln7_2, = ax.plot((pos_botom[0], pos_l[0]), (-pos_botom[1], -pos_l[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                ln8_2, = ax.plot((pos_top[0], pos_r[0]), (-pos_top[1], -pos_r[1]), linestyle=linestyle, color=linecolor, linewidth=linewidth, solid_capstyle='butt')
                
                if symmetric == False:
                    if map_order == 0:
                        ln7.remove()
                        ln8.remove()
                        ln7_2.remove()
                        ln8_2.remove()
                    if map_order == 1:
                        ln5.remove()
                        ln6.remove()
                        ln5_2.remove()
                        ln6_2.remove()
                        
def _mark_contact_regions(row_regions, 
                         markRegions, 
                         binsize,
                         ax, 
                         map_type, 
                         map_order,
                         symmetric, 
                         linestyle, 
                         linecolor, 
                         linewidth,
                         only_cis):
    
    plot_region_df = _get_pets(row_regions, binsize)
 
    mark_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in markRegions])
    mark_GenomeRegions_cp = mark_GenomeRegions.copy()
    for i, row in mark_GenomeRegions_cp.iterrows():
        plot_region_df[i] = 0
        bin_s = int(row['fetch_start']/binsize)
        if row['fetch_start']%binsize > 0:
            bin_s = bin_s + 1

        bin_e = int(row['fetch_end']/binsize)
        if row['fetch_end']%binsize > 0:
            bin_e = bin_e + 1
        one_region = plot_region_df[plot_region_df["chrom"]==row['chrom']]
        one_region = one_region[(one_region['chrbin']>=bin_s) & (one_region['chrbin']<=bin_e)]
        if one_region.shape[0] > 0:
            plot_region_df.loc[one_region.index, i] = 1

    if map_type == "square":
        _mark_contact_regions_square(ax, markRegions, plot_region_df, map_order, symmetric, linewidth, linestyle, linecolor, only_cis)
    elif map_type == "triangle":
        _mark_contact_regions_triangle(ax, markRegions, plot_region_df, map_order, symmetric, linewidth, linestyle, linecolor, only_cis)
    else:
        pass
       

def mapc_markline(ax: Optional[Axes] = None,
                  row_regions: Union[pd.DataFrame, None] = None,
                  mark_regions: Union[Sequence[str], str, None] = None,
                  binsize: Union[int, None] = 10000,
                  map_order: int = 0,
                  map_type: Union[str, None] = 'squ',
                  symmetric: bool = False,
                  trans_ax: bool = False,
                  only_cis: bool = False,
                  linestyle: Union[str, None]='--',
                  linecolor: Union[str, None]='k',
                  linewidth: Union[int, None]=1,
                  show_regions_edge: bool = False,
              ):
    """
    The ``mapc_markline`` function used for mark contact regions on the heatmap, 
    It should be used in conjunction with ``mapc``. The parameters "symmetric" and 
    "trans_ax" which should be consistent with the parameters used in ``mapc``.
    
    Parameters
    ----------
    ax: :class:`matplotlib.axes.Axes` object
        utilizes the same `Axes` as used in ``mapc``
    row_regions: `pd.DataFrame`
        The `row_regions` attribute of the object returned by ``tc.tl.extractContactRegions``
        represents the regions extracted from the rows.
    mark_regions: `np.ndarray`
        matrix for plot upper or right of heatmap
    binsize: `int`
        binsize of raw heatmap
    map_order: `int`
        one option of [0, 1]
        The `map_order` parameter in the context of `trackc.mapc` is associated with the `mat` and `mat2` parameters.
        If use `trackc.mapc` to generate a heatmap in the top-right quadrant, then for "map_order" here, you should set it to 0.
        If you use "mapc" to generate a heatmap in the bottom-left quadrant, then for "map_order" here, you should also set it to 0.
    map_type: `str`
        same as the `trackc.mapc` map_type parameter.
        optional is ['square', 'squ', 'triangle', 'tri', 'rectangle', 'rec'], default is `square`
        que is the same as square, tri is the same as triangle, rec is the same as rectangle.
    symmetric: `bool`
        same as the `trackc.mapc` symmetric parameter.
    trans_ax: `bool`
        same as the `trackc.mapc` trans_ax parameter.
    only_cis: `bool`
        If set to True, the trans regions will not be displayed.
    linestyle: `str`
        mark line style
    linecolor: `str`
        mark line color
    linewidth: `str`
        mark line width
    show_regions_edge: `bool`
        Whether to display the spliced region as input from the `row_regions` parameter

    """
    if map_order == 1:
        for child in ax.get_children():
            if isinstance(child, Axes):
                ax = child
                break
        

    if map_type=="squ":
        map_type="square"
    if map_type=="tri":
        map_type="triangle"
    if map_type=="rec":
        map_type="rectangle"

    if mark_regions != None:
        _mark_contact_regions(row_regions, mark_regions, binsize, ax, map_type, map_order, symmetric, linestyle, linecolor, linewidth, only_cis)

    if show_regions_edge==True:
        _plot_contact_regions_line(ax, row_regions, map_order, map_type, trans_ax, linewidth, linestyle, linecolor)

def _plot_contact_regions_line(ax, 
                               rl_regions, 
                               map_order=0,
                               map_type='triangle',
                               trans_ax=False,
                               linewidth=1,
                               linestyle='--',
                               linecolor='k'):
    if trans_ax == True:
        print('trans_ax should be False')

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

        if map_type == 'triangle':
            if map_order==0:
                xA, yA, xB, yB = v_sum/2, v_sum, v_sum, 0
                xA2, yA2, xB2, yB2 = v_sum, 0, v_sum+(full_len-v_sum)/2, full_len-v_sum
            else:
                xA, yA, xB, yB = v_sum/2, -v_sum, v_sum, 0
                xA2, yA2, xB2, yB2 = v_sum, 0, v_sum+(full_len-v_sum)/2, -(full_len-v_sum)
        
        if map_type == 'square':
            if map_order==0:
                xA, yA, xB, yB = v_sum-0.5, -0.5, v_sum-0.5, v_sum-0.5
                xA2, yA2, xB2, yB2 = v_sum-0.5, v_sum-0.5, full_len-0.5, v_sum-0.5
            else:
                xA, yA, xB, yB = -0.5, v_sum-0.5, v_sum-0.5, v_sum-0.5
                xA2, yA2, xB2, yB2 = v_sum-0.5, v_sum-0.5, v_sum-0.5, full_len-0.5

        ax.plot([xA, xB], [yA, yB], linestyle, linewidth=linewidth, color=linecolor, solid_capstyle='butt')
        ax.plot([xA2, xB2], [yA2, yB2], linestyle, linewidth=linewidth, color=linecolor, solid_capstyle='butt')



