#regions := chr8:127000000-129200000 chr14:96500000-99300000
r1=chr8:127000000-129200000
r2=chr14:96500000-99300000
#regions := chr14:96500000-99300000
binsize=10000
chrom_sizes=/Users/yuanzan/Documents/ref/chrom_GRCh38.105.sizes
in_cool=GSM4604287_1360.iced.mcool::/resolutions/10000
out_cool=/Users/yuanzan/Documents/github/seqyuan/trackc_data/tutorials/4C/GSM4604287_1360.sub.cool

subset:
	@if [ -e cool.tmp.txt ]; then \
		rm cool.tmp.txt;\
	fi
	@$(foreach regin, $(regions), cooler dump --join -r $(regin) $(in_cool) >>cool.tmp.txt;)
	cat cool.tmp.txt |cooler load --format bg2 $(chrom_sizes):$(binsize) - $(out_cool)
	#rm cool.tmp.txt

extract:
	@if [ -e cool.tmp.txt ]; then \
		rm cool.tmp.txt;\
	fi
	cooler dump --join -r $(r1) $(in_cool) >cool.tmp.txt
	cooler dump --join -r $(r2) $(in_cool) >>cool.tmp.txt
	cooler dump --join -r $(r1) -r2 $(r2) $(in_cool) >>cool.tmp.txt
	cat cool.tmp.txt |cooler load --format bg2 $(chrom_sizes):$(binsize) - $(out_cool)
	rm cool.tmp.txt

