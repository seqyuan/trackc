from trackc.tl._getRegionsCmat import extractContactRegions
from .mapc import mapC, getData2Map
from matplotlib.axes import Axes
from typing import Union, Sequence, Optional
from matplotlib.colors import Colormap
import pandas as pd
import numpy as np
import cooler
import seaborn as sns
fruitpunch = sns.blend_palette(['white', 'red'], as_cmap=True)

def _plot4C_line_bar(data,
                     ax,
                     color,
                     track_type,
                     logdata,
                     trim_range,
                     minrange,
                     maxrange,
                     ):
    data, maxrange, minrange = getData2Map(data, maxrange=maxrange, minrange=minrange, trim_range=trim_range, inplace=False)
    print("maxrange:", maxrange ,"minrange:",minrange)
    
    if logdata:
        data = np.log2(data)
    
    if track_type=='line':
        ax.plot(range(len(data[0,:])), data[0,:], color=color, solid_capstyle='butt')

    if track_type=='bar':
        ax.bar(x=range(len(data[0,:])), height=data[0,:], bottom=minrange, width=1, color=color, align='edge')
    return (data, maxrange, minrange)

def _get_pets(df, binsize=10000):
    pets = df.reset_index()
    pets_df = pd.DataFrame(index=range(pets['cbins'].sum()))
    pets_df['chrom'] = 'None'
    pets_df['bin_N'] = 0
    
    chroms = []
    bins = []
    
    for i, row in pets.iterrows():
        chroms.extend([row['chrom']] * row['cbins'])
         
        bin_s = int(row['fetch_start']/binsize)
        #if row['fetch_start']%binsize > 0:
        #    bin_s = bin_s

        bin_e = int(row['fetch_end']/binsize)
        if row['fetch_end']%binsize > 0:
            bin_e = bin_e + 1
        
        chrbins = list(range(bin_s, bin_e))
        
        if row['isReverse'] == True:
            chrbins = chrbins[::-1]
        bins.extend(chrbins)

    pets_df['chrom'] = chroms
    pets_df['chrbin'] = bins
    return pets_df


def virtual4C(ax: Optional[Axes] = None,
              clr: cooler.Cooler = None,
              balance: bool = False,
              #divisive_weights = None,
              target: Union[str, None] = None,
              contact_regions: Union[Sequence[str], str, None] = None, 
              track_type: Union[str, None] = 'line',
              color: Union[str, None] = 'tab:blue',
              cmap: Union[Sequence[Colormap], str, None] = fruitpunch,
              target_color: Union[str, None] = None,
              target_name: Union[str, None] = None,
              logdata: Union[Sequence[bool], bool] = False, 
              trim_range: float = 0.98,
              minrange: float = None,
              maxrange: float = None,
              label: Optional[str] = None,
              label_rotation: Union[int, None] = 0,
              label_fontsize: Optional[int] = 12,
              ):
    data = extractContactRegions(clr, balance=balance, row_regions=target, col_regions=contact_regions)
    
    if track_type in ['line', 'bar']:
        _, maxrange, minrange = _plot4C_line_bar(data=data.cmat,
                ax=ax,
                color=color,
                track_type = track_type,
                logdata=logdata, 
                trim_range=trim_range, 
                minrange=minrange, 
                maxrange=maxrange, 
                )
        
        cols = _get_pets(data.col_regions, clr.binsize)
        rows = _get_pets(data.row_regions, clr.binsize)
        cols['target'] = 0
        cols.loc[cols['chrbin'].isin(rows['chrbin']), 'target'] = 1
        targets_point = cols.query('target==1')
        target_bar = ax.bar(x=targets_point.index[0], height=maxrange, 
               bottom=minrange, width=1, color=target_color, align='edge', 
               label=target_name)
        #ax.bar_label(target_bar, label_type='edge')
        ax.text(targets_point.index[0], minrange+(maxrange-minrange)/2, target_name, va='center', ha='right', rotation=90)
        ax.set_ylim(minrange, maxrange)

    if track_type == 'heatmap':
        mapC(mat=data.cmat, cmap=cmap, logdata=logdata, 
             trim_range=trim_range, minrange=minrange, maxrange=maxrange, 
             ax=ax, map_type='square', ax_on=False)
        ax.set_aspect(aspect='auto')
    
    
    

    ax.tick_params(bottom =False,top=False,left=False,right=True)
    ax.yaxis.tick_right()
    #ax.set_xticklabels("")
    ax.set_xticklabels("")
    ax.set_xlim(0, data.cmat.shape[1]-1)
    ax.set_ylabel(label, fontsize=label_fontsize, rotation=label_rotation, ha='right', va='center')




