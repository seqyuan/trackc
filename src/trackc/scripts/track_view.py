import pandas as pd
import numpy as np
import os
import sys,re

import matplotlib
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
matplotlib.use('Agg')
#%matplotlib inline

import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pyBigWig
import yaml

import warnings
warnings.filterwarnings('ignore')


def plot_bwTrack(ax, bw, ylabel, chrom, start, end, resolution=20000 , yminx=5, ymaxx=95,rotation=0, fl='%0.2f',color='#464451'):
    ax.tick_params(bottom =False,top=False,left=True,right=False)
    ax.spines['left'].set_color('k')
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_color('k')
    ax.spines['bottom'].set_linewidth(0.5)

    x = int((end-start)/resolution)

    plot_list = bw.stats(chrom, start, end, type="mean", nBins=x)
    plot_list = [0 if v is None else v  for v in plot_list]
    
    width = 1
    ax.bar(x=range(0,x), height=plot_list, width=1, bottom=[0]*(x),color=color,align="edge",edgecolor=color)    
    ax.set_xlim(0,x)

    #ax.set_xticks([])
    #ax.set_xticklabels([])

    ymin = np.percentile(plot_list,yminx)
    ymax = np.percentile(plot_list,ymaxx)

    ax.set_yticks([ymin, ymax])
    ax.set_yticklabels([fl % ymin, fl % ymax], fontsize=7)
    
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(ylabel, fontsize=8, rotation=rotation, horizontalalignment='right',verticalalignment='center')
    
    ax.set_xticks([])
    ax.set_xticklabels('')
    
    
def plot_bwTracks(tracks, outfile_prefix="track_test"):
    fig = plt.figure(figsize=(int(tracks["figure"]["size"]["width"]), int(tracks["figure"]["size"]["height"])), facecolor='white')
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.8, hspace=0.2, wspace=0.15)

    gs = fig.add_gridspec(int(tracks["figure"]["row"]), int(tracks["figure"]["col"]))
    ax = fig.add_subplot(gs[0, 0], facecolor='white')
    
    for i, track in enumerate(tracks["tracks"]):
        if track["gs"]["row"]["to"] == None:
            if track["gs"]["col"]["to"] == None:
                ax = fig.add_subplot(gs[int(track["gs"]["row"]["from"]), int(track["gs"]["col"]["from"])], facecolor='white')
            else:
                ax = fig.add_subplot(gs[int(track["gs"]["row"]["from"]), int(track["gs"]["col"]["from"]):int(track["gs"]["col"]["to"])], facecolor='white')
        else:
            if track["gs"]["col"]["to"] == None:
                ax = fig.add_subplot(gs[int(track["gs"]["row"]["from"]):int(track["gs"]["row"]["to"]), int(track["gs"]["col"]["from"])], facecolor='white')
            else:
                ax = fig.add_subplot(gs[int(track["gs"]["row"]["from"]):int(track["gs"]["row"]["to"]), int(track["gs"]["col"]["from"]):int(track["gs"]["col"]["to"])], facecolor='white')
        
        if track["type"] == "bw":
            plot_bwTrack(ax, track["bwData"], track["name"], track["pars"]["chrom"], track["pars"]["start"], track["pars"]["end"], resolution=track["pars"]["resolution"], color=track["pars"]["color"])                                                                                                    

    #plt.show()
    fig.savefig(outfile_prefix + ".pdf")
    fig.savefig(outfile_prefix + ".png")
    

def get_yaml_data(yaml_file):
    file = open(yaml_file, 'r', encoding="utf-8")
    file_data = file.read()
    file.close()

    data = yaml.load(file_data)
    return data

def usage():
    print("python3 track_view.py config.yml outdir outfile_prefix")
    sys.exit(1)

def main(args):
    yml_file, outdir, outfile_prefix =  args
    tracks = get_yaml_data(yml_file)

    for i, track in enumerate(tracks["tracks"]):
        if track["type"] == "bw":
            tracks["tracks"][i]["bwData"] = pyBigWig.open(track["pars"]["bwfile"])

    outfile_prefix = os.path.join(outdir,outfile_prefix)
    print(outfile_prefix)
    plot_bwTracks(tracks, outfile_prefix)



if __name__ == '__main__':
    if len(sys.argv) != 4:
        usage()
    main(sys.argv[1:])
