import sys
#sys.path.append("/Users/yuanzan/Documents/github/seqyuan/trackC")
from trackc import heatmap

import matplotlib.pyplot as plt

import matplotlib
#%matplotlib inline
matplotlib.use('TkAgg')


import yaml
import os

def get_yaml_data(yaml_file):
    file = open(yaml_file, 'r', encoding="utf-8")
    file_data = file.read()
    file.close()

    data = yaml.load(file_data)
    return data

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

def plot_track_map(tracks, outfile="heatmap.pdf"):
    #fig = plt.figure(figsize=(defaults["figure"]["size"]["width"], defaults["figure"]["size"]["height"]), constrained_layout=False, facecolor='white')
    fig = plt.figure(figsize=(10, 12), constrained_layout=False, facecolor='white')
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.8, top=0.8, hspace=0.1, wspace=0.15)

    gs = fig.add_gridspec(int(tracks["figure"]["row"]), int(tracks["figure"]["col"]))
    ax = gs[0, 0]
    
    maxr = None
    minr = None
    
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
        
        im, pm, pp = heatmap.triViewC(ax, tracks["tracks"][i]["plotMat"], ylabel=tracks["tracks"][i]["name"], cmap=heatmap.colorC(cname="RdBu_r",bottom_color="#09188e"),trim_range=0.98)                                                                                                    

        if i == 0:
            heatmap.colorbar_upright(ax,im,pp,pm)                                             
    ax0 = fig.add_subplot(gs[0:2, 0], facecolor='None')                                                                                                        
    plot_xticks(ax0, "chr7", 17000000, 28000000)

    fig.savefig(outfile)   


yaml_path = "config.yml"
tracks = get_yaml_data(yaml_path)

for i, track in enumerate(tracks["tracks"]):


    if track["type"] == "triView":
        mat = heatmap.read_chrom_matrix(track["pars"]["chrom"], track["pars"]["abs.bed"], track["pars"]["matrix"])
        tracks["tracks"][i]["TriangularMat"] = heatmap.mat2Triangular(mat)
        plot_mat, s, e = heatmap.subset2triMat(tracks["tracks"][i]["TriangularMat"], start=track["pars"]["start"], end=track["pars"]["end"], resolution=track["pars"]["resolution"], show_distance=track["pars"]["show_distance"])
        tracks["tracks"][i]["plotMat"] = plot_mat
        
plot_track_map(tracks)   



