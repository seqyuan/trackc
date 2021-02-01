import numpy as np
import pandas as pd

import matplotlib
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
plt.style.use('ggplot')
from scipy import sparse
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import warnings
warnings.filterwarnings('ignore')


def read_chrom_matrix(chrom, bed_file, mat_file):
    #load HiC-Pro resolution_abs.bed and matrix return one chrom N x N matrix
    bed = pd.read_table(bed_file, header=None,index_col = None,names=['chrom','start','end','bNum'],encoding='utf-8')
    bed = bed[bed['chrom'] == chrom]
    minN = bed.iloc[0,3]
    maxN = bed.iloc[-1,3]
    N = maxN - minN + 1

    df = pd.read_table(mat_file, header=None,index_col = None,names=['b1','b2','contact'],encoding='utf-8')
    df= df[(df['b1'] >= minN) & (df['b1'] <= maxN) & (df['b2'] >= minN) & (df['b2'] <= maxN)]
    df['b1'] = df['b1'] - minN
    df['b2'] = df['b2'] - minN
    bed['bNum'] = bed['bNum'] - minN

    counts = sparse.coo_matrix((df['contact'], (df['b1'], df['b2'])), shape=(N, N), dtype=float)
    df = pd.DataFrame(counts.todense(),index = bed['bNum'].tolist(),columns = bed['bNum'].tolist())
    return df

def mat2Triangular(raw_mat=None):
    # convert N x N matrix to Triangular for plot heatmap
    mat=np.array(raw_mat)
    n=0;length=len(mat)
    tri_mat=np.zeros([length,length*2])
    tri_mat[tri_mat==0]=np.nan
    for i in range(length):
        curl=np.array(np.diag(mat,i))
        tri_mat[i,n:(length*2-n)]=curl.repeat(2)
        n+=1
    return tri_mat

def subset2triMat(mat,start=17000000,end=28940000,resolution=40000,show_distance=2500000):
    distance = int(show_distance/resolution)
    start_bin = int(start/resolution)
    end_bin = int(end/resolution)
    if end%resolution > 0:
        end_bin = end_bin + 1
        
    plot_mat = mat[0:distance, start_bin*2:end_bin*2]
    return plot_mat, start_bin*2, end_bin*2

    
def getData2Map(mat, maxrange=None, minrange=None, logdata=True, trim_range=98):
    if logdata:
        mat = np.log2(mat)
        
    mat[mat == np.inf] = 0
    mat[mat == -np.inf] = 0    
        
    df = pd.DataFrame(np.ravel(mat))
    #df = df[df != np.inf]
    #df = df[df != -np.inf]
    df = df[df != 0]
    
    if np.nanmax(mat)>0:
        if trim_range <100 and maxrange==None and minrange==None:
            xmaxrange = np.nanpercentile(abs(df), 100-trim_range)
            xminrange = np.nanpercentile(abs(df), trim_range)
        else:
            if maxrange==None:
                xmaxrange = abs(df).max()[0]
            else:
                xmaxrange = maxrange
            if minrange==None:
                xminrange = abs(df).min()[0]
            else:
                xminrange = minrange

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

def colorC(cname="RdBu_r", bottom_color="#ffffff", bad_color="white",over_color="white",under_color="white",alpha=0):
    #cname = ["twilight_shifted", "jet", "RdBu_r", "RdGy_r", "BrBG_r", "hot_r", "Spectral_r"]
    #cname = ["terrain_r", "ocean", "gist_earth", "gist_stern_r", "tab20b", "twilight"]
    from matplotlib.colors import ListedColormap
    cmap=plt.get_cmap(cname)

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

def colorC(cname="RdBu_r", bottom_color='#ffffff', bad_color="white",over_color="white",under_color="white",alpha=0):
    #cname = ["twilight_shifted", "jet", "RdBu_r", "RdGy_r", "BrBG_r", "hot_r", "Spectral_r"]
    #cname = ["terrain_r", "ocean", "gist_earth", "gist_stern_r", "tab20b", "twilight"]
    from matplotlib.colors import ListedColormap
    cmap=plt.get_cmap(cname)

    cmap.set_bad(color=bad_color, alpha=alpha)
    cmap.set_over(color=over_color, alpha=alpha)
    cmap.set_under(color=under_color, alpha=alpha)

    if bottom_color==None:
        return cmap
    else:
        bottom_color = hex2rgb(bottom_color)
    
        newcolors = cmap(np.linspace(0, 1, 256))
        #white = np.array([1, 1, 1, 1])
        newcolors[0, :] = bottom_color
        newcmap = ListedColormap(newcolors)
        return newcmap


def triViewC(ax, mat, ylabel='sample',aspect='auto',cmap=plt.get_cmap("RdBu_r"), logdata=True, trim_range=0.98,maxrange=None,minrange=None):
    #aspect=['auto','equal']
    mat, maxrange, minrange = getData2Map(mat, maxrange=maxrange, minrange=minrange, logdata=logdata, trim_range=trim_range)
    print("maxrange:", maxrange ,"minrange:",minrange)
    im = ax.imshow(mat, cmap=cmap,interpolation="nearest",origin='lower', clim=(minrange,maxrange))
    ax.tick_params(bottom =False,top=False,left=False,right=False) #去掉tick线
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(ylabel)
    return im, minrange, maxrange

def colorbar_upright(axm,im, vmax, vmin):
    axins1 = inset_axes(axm, width="8%", height="12%", loc=2, bbox_to_anchor=(0.8, 0.3, 2.3, 1), bbox_transform=axm.transAxes)
    ticks = [vmin, vmax]
    cbar=plt.colorbar(im, cax=axins1, orientation='horizontal', ticks=ticks)
    cbar.ax.tick_params(axis='x', pad=2)
    cbar.ax.set_xticklabels(['Low','High'],rotation=0,fontsize=10)
    cbar.ax.xaxis.tick_top()
    cbar.ax.tick_params(bottom =False,top=False,left=False,right=False)
    
def colorbar_right(axm,im, vmax, vmin):
    axins1 = inset_axes(axm, width="0.7%", height="80%", loc=3, bbox_to_anchor=(1, 0.2, 2.3, 0.8), bbox_transform=axm.transAxes,)
    ticks = [vmin, vmax]
    cbar=plt.colorbar(im, cax=axins1, orientation='vertical', ticks=ticks)
    cbar.ax.set_yticklabels(['Low','High'],rotation=0,fontsize=10)
    cbar.ax.tick_params(bottom =False,top=False,left=False,right=False)


def set_ticks(start, end):
    # for plot_xticks
    distance = end - start
    if distance >= 2000000:
        ticks = list(range(int(start/1000000), int(end/1000000)))[1:]
        tick_labels = [str(x) for x in ticks]
        tick_labels[-1] = tick_labels[-1] + 'Mb'
        ticks = [x*1000000 for x in ticks]
        return ticks, tick_labels
    
    if distance >= 500000:
        ticks = list(range(int(start/1000), int(end/1000)))[1:]
        tick_labels = [str(x) for x in ticks]
        tick_labels[-1] = tick_labels[-1] + 'Kb'
        ticks = [x*1000 for x in ticks]
        return ticks, tick_labels

def plot_xticks(ax, chrom, start, end):
    ax_ticks, ax_tick_labels = set_ticks(start, end)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    for tick,label in zip(ax_ticks, ax_tick_labels):
        ax.text(tick, 1.5, label, fontsize=10, ha='center', va='top')
        ax.plot([tick,tick],[2.5,3.5],'k-',lw=0.5)
        
    ax.axhline(y=2.5, c="k", ls="-", lw=1.5)
    ax.text(start, 3.5,"{0}:{1}-{2}".format(chrom, start, end), fontsize=15)
    ax.set_xlim([start, end])
    ax.set_ylim([-10,10])

def plot_bottom_label(ax, plot_start, plot_end, tick_start, tick_end, gene_start, gene_end, gene_name):
    ax.plot([tick_start, tick_start],[0,2],'r-',lw=0.5)
    ax.plot([tick_end, tick_end],[0,2],'r-',lw=0.5)
    
    ax.set_xlim([plot_start, plot_end])
    ax.set_ylim([10,0])

    ax.bar(gene_start,3,width=gene_end-gene_start,bottom=[0],color='red',align='edge')
    #ax.bar([gene_start, gene_end],[0,0],'r-',lw=10)
    name_x = (gene_start + gene_end)/2
    ax.text(name_x, 1, gene_name, fontsize=10, ha='center', va='top')

    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])

