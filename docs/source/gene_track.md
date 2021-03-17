# gene_track
## Description
A track for gene bed12 files.

## Parameters
### Necessary:
* gene bed12 file

gene_bed12 colnames: 
["chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]

follow colnames are necessary to plot gene track: "chrom","start","end","name" "blockSizes","blockStarts"

### Optional:
* ax: matplotlib.axis.Axis
* gene_bed: pandas.DataFrame, colnames = ["chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
* chrom: str, chromosome name, example: chr1
* start: int, plot_start point, example: 30000000
* end: int, plot_end point, example: 35000000
* line: int, 1 (default), to avoid genetrack overlapping, plot genetrack lines, 
* gene_col: list, \[\] (default), genes mark a special color(red), example: ["DACH1", "PIBF1"]
* gene_name_fontszie: int, 5 (default), gene name label fontszie
