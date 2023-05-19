from trackc.tl._getRegionsCmat import GenomeRegion
from matplotlib.axes import Axes
from typing import Union, Sequence, Optional
import pandas as pd
from trackc.pa import trackcl_11

def zoomin(ax: Optional[Axes] = None,
           row_regions: Union[pd.DataFrame, None] = None,
           zoomin_regions: Union[Sequence[str], str, None] = None,
           colors: Union[Sequence[str], None] = None,
           alpha: float = 1,
           line_on: bool = False,
           fill: bool = True,
           ):
    """\
    map_order:
        mapc mat or mat2
    """
    
    row_GRs = _region_pos(row_regions)
    zoomin_GRs = _region_pos(zoomin_regions)
    zoomin_GRs['ors'] = None
    zoomin_GRs['ore'] = None
    zoomin_GRs_cp = zoomin_GRs.copy()
    for i, row in zoomin_GRs_cp.iterrows():
        start, end = _get_zoom_origin_pos(row['chrom'], row['fetch_start'], row['fetch_end'], row_GRs)
        ss = start
        ee = end
        if row['isReverse'] == True:
            ss = end
            ee = start
        zoomin_GRs.loc[i, 'ors'] = ss
        zoomin_GRs.loc[i, 'ore'] = ee
    
    zoomin_GRs = zoomin_GRs.reset_index()

    if colors == None:
        colors = trackcl_11
    
    for i, row in zoomin_GRs.iterrows():
        x = row[['ps', 'pe', 'ore', 'ors']]
        y = [0,0,1,1]

        if fill == True:
            ax.fill(x, y, color=colors[i], alpha=alpha)

        if line_on==True:
            ax.plot((x[0], x[3]), (y[0], y[3]),color=colors[i], alpha=alpha, solid_capstyle='butt')
            ax.plot((x[1], x[2]), (y[1], y[2]),color=colors[i], alpha=alpha, solid_capstyle='butt')
        
    ax.set_axis_off()
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])

def _region_pos(regions):
    GRs = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    GRs['len'] = GRs['fetch_end'] - GRs['fetch_start']
    GRs['pos_e'] = GRs['len'].cumsum()
    GRs['pos_s'] = GRs['pos_e'] - GRs['len']
    GRs_fulllen = GRs['len'].sum()
    GRs['ps'] = GRs['pos_s']/GRs_fulllen
    GRs['pe'] = GRs['pos_e']/GRs_fulllen
    return GRs


def _get_zoom_origin_pos(zoom_chrom, zooms, zoome, origin_pos_df):
    fulllen = origin_pos_df['len'].sum()
    origin_pos = origin_pos_df[origin_pos_df['chrom']==zoom_chrom]
    origin_pos = origin_pos[(zooms>=origin_pos['fetch_start']) & (zoome<=origin_pos['fetch_end'])]
    
    start = None
    end = None
    if origin_pos.shape[0] == 0:
        print('{0}:{1}-{2} not in origion regions'.format(zoom_chrom, zooms, zoome))
    else:
        ix = origin_pos.index[0]
        start_dis = origin_pos_df.loc[ix, 'pos_s'] + abs(zooms - origin_pos_df.loc[ix, 'start'])
        end_dis = origin_pos_df.loc[ix, 'pos_s'] + abs(zoome - origin_pos_df.loc[ix, 'start'])
        start = start_dis/fulllen
        end = end_dis/fulllen
    return start, end


