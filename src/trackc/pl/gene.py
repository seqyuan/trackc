from matplotlib.axes import Axes
from typing import Union, Optional, Sequence
import pandas as pd
from trackc.tl._getRegionsCmat import GenomeRegion
from .bigwig import make_multi_region_ax


def gene_track(gene_bed: pd.DataFrame,
               ax: Optional[Axes] = None,
               regions: Union[Sequence[str], str, None] = None,
               track_type: Union[str, None] = 'gene',
               pos_strand_gene_color: Union[str, None] = '#3366CC',
               neg_strand_gene_color: Union[str, None] = '#EECFA1',
               line: Union[int, None] = 1,
               gene_fontszie: Union[int, None] = 5,
               ylabel: Optional[str] = None,
               label_rotation: Union[int, None] = 0,
               label_fontsize: Optional[int] = 12,
               ):
    """\
    Plot multi-regions gene track.
    
    Parameters
    ----------
    gene_bed
        ``pd.DataFrame``:
    ax

    regions

    track_type

    pos_strand_gene_color

    neg_strand_gene_color

    line

    gene_fontszie

    ylabel

    label_rotation

    label_fontsize

    """

    if isinstance(regions, list):
        line_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in regions])
    else:
        line_GenomeRegions = GenomeRegion(regions).GenomeRegion2df()

    axs = make_multi_region_ax(ax, line_GenomeRegions)
    line_GenomeRegions = line_GenomeRegions.reset_index()

    ax.set_ylabel(ylabel, fontsize=label_fontsize, rotation=label_rotation, horizontalalignment='right',verticalalignment='center')
    spines = ['top', 'bottom', 'left', 'right']
    for i in spines:
        ax.spines[i].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels('')
    ax.set_yticks([])
    ax.set_yticklabels('')  

    gene_bed.columns =['chrom','start','end','name',"score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
    for ix, row in line_GenomeRegions.iterrows(): 
        if track_type == "gene":
             plot_gene(axs[ix], gene_bed, row['chrom'], row['fetch_start'], row['fetch_end'], needReverse=row['isReverse'], pos_strand_gene_color=pos_strand_gene_color, neg_strand_gene_color=neg_strand_gene_color, line=line, fontszie=gene_fontszie)
    
def plot_gene(ax, gene_bed, chrom, start, end, needReverse=False, pos_strand_gene_color='#3366CC', neg_strand_gene_color='#EECFA1', line=1, fontszie=5):
    gene_bed = gene_bed[gene_bed['chrom']==chrom]
    gene_bed_plot = gene_bed[((gene_bed['start'] >= start) & (gene_bed['start'] <= end)) | ((gene_bed['end'] >= start) & (gene_bed['end'] <= end))]
    gene_bed_plot = gene_bed_plot.sort_values(by='end')
    #print(gene_bed_plot
    
    plot_gene_num = gene_bed_plot.shape[0]

    ii = 0
    for i,row in gene_bed_plot.iterrows():
        #col = pos_strand_gene_color
        text_col = pos_strand_gene_color
        
        if row["strand"] == "-":
            #col = neg_strand_gene_color
            text_col = neg_strand_gene_color
        
        #text_col = col
        plot_y = ii%line
        
        ax.plot((row['start'], row['end']), (plot_y + 0.5, plot_y+0.5), color='k', linewidth=1, solid_capstyle='butt')
        starts = [int(x) for x in row["blockStarts"].split(",")]
        widths = [int(x) for x in row["blockSizes"].split(",")]
        
        ax.bar(x=starts, height=0.4, width=widths, bottom=plot_y+0.3, \
                edgecolor='k', linewidth=1, align='edge', color='k')
        
        if row['start'] < start:
                row['start'] = start
        if row['end'] > end:
                row['end'] = end
        
        arrow_s = row['end']
        dx = 0.1
        if row["strand"] == "-":
            arrow_s = row['start']
            dx = -0.1
        ax.arrow(arrow_s, plot_y+0.5, dx, 0, 
            overhang=0, width=0,
            head_width=0.25,
            head_length=10000,
            length_includes_head=False,
            color=text_col,
            linewidth=1)

        if (row['name'] in gene_bed_plot.iloc[-4:,:]['name']) or (int(line/(ii+1)) < 2):
            ha = 'right'
            genename = row['name'] + "  "
            xpos = row['start']
            if needReverse:
                ha = 'left'
                genename = "  " + row['name']
                #xpos = row['end']
            ypos = plot_y + 0.5
            if line == 1:
                 xpos = row['start'] + abs(row['start']-row['end'])/2
                 ypos = 0.8
                 ha = 'center'
            ax.text(xpos, ypos, genename + "  ", ha=ha, va='center',color=text_col, fontsize=fontszie)
        else:
            ha = 'left'
            genename = "  " + row['name']
            xpos = row['end']
            if needReverse:
                ha = 'right'
                genename = row['name'] + "  "
                #xpos = row['start']

            ypos = plot_y + 0.5
            if line == 1:
                 xpos = row['start'] + abs(row['start']-row['end'])/2
                 ypos = 0.8
                 ha = 'center'
            ax.text(xpos, ypos, genename, ha=ha, va='center',color=text_col, fontsize=fontszie)

        ii+=1
            
    xlim_s = start
    xlim_e = end
    if needReverse == True:
         xlim_s = end
         xlim_e = start

    ax.set_xlim(xlim_s, xlim_e)
    ax.set_ylim(top=0, bottom=line)
    if plot_gene_num < line:
        ax.spines['bottom'].set_position(('data', plot_gene_num))
    
    for i in ['left','top','right']:
        ax.spines[i].set_color('none')
        ax.spines[i].set_linewidth(0)
    ax.spines["bottom"].set_color('black')
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(bottom =True,top=False,left=False,right=False)
    ax.set_xticklabels("")
    ax.set_yticklabels("")

