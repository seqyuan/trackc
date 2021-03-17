def gene_track(ax, gene_bed, chrom, start, end, line=1, gene_col=[], gene_name_fontszie=5):
    # gene_bed colnames["chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
    # follow colnames are necessary
    # "chrom","start","end","name" "blockSizes","blockStarts"
    gene_bed = gene_bed[gene_bed['chrom']==chrom]
    gene_bed_plot = gene_bed[((gene_bed['start'] >= start) & (gene_bed['start'] <= end)) | ((gene_bed['end'] >= start) & (gene_bed['end'] <= end))]
    gene_bed_plot = gene_bed_plot.sort_values(by='end')
    #print(gene_bed_plot
    
    plot_gene_num = gene_bed_plot.shape[0]
    print(plot_gene_num)
    if line == 1:
        for i,row in gene_bed_plot.iterrows():
            ax.plot((row['start'], row['end']), (0.5, 0.5), color='b')
            starts = [int(x) for x in row["blockStarts"].split(",")]
            widths = [int(x) for x in row["blockSizes"].split(",")]
            
            ax.bar(x=starts, height=1, width=widths, bottom=0, \
                   edgecolor='black', linewidth=0, align='edge', color='b', ecolor=None)
    else:
        ii = 0
        for i,row in gene_bed_plot.iterrows():
            col = "#E69F00"
            if row["strand"] == "-":
                col = "#56B4E9"
            
            #text_col = col
            text_col = 'k'
            if row["name"].isin(gene_col):
                text_col = "red"
                
            plot_y = ii%line
            ax.plot((row['start'], row['end']), (plot_y + 0.5, plot_y+0.5), color=col)
            starts = [int(x) for x in row["blockStarts"].split(",")]
            widths = [int(x) for x in row["blockSizes"].split(",")]
            
            ax.bar(x=starts, height=1, width=widths, bottom=plot_y, \
                   edgecolor=col, linewidth=0.3, align='edge', color=col)
            
            if (gene_bed_plot.iloc[-1]['name'] == row['name']) or (gene_bed_plot.iloc[-2]['name'] == row['name']) or (gene_bed_plot.iloc[-3]['name'] == row['name']):
                ax.text(row['start'] - 5000, plot_y + 0.5, row['name'], ha='right', va='center',color=text_col, fontsize=fontszie)
            else:
                ax.text(row['end'] + 5000, plot_y + 0.5, row['name'], ha='left', va='center',color=text_col, fontsize=fontszie)

            ii+=1
    ax.set_xlim(start,end)
    for i in ['bottom','left','top','right']:
        ax.spines[i].set_color('black')
        ax.spines[i].set_linewidth(0.5)
    ax.tick_params(bottom =True,top=False,left=False,right=False)
    #ax.set_xticklabels("")
    ax.set_yticklabels("")

    