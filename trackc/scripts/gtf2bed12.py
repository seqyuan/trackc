class _genebed12(object):
    """bed12 object"""

    def __init__(self, chrom, start, end, name, strand, gene_biotype):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = 0
        self.strand = strand
        self.thickStart = end
        self.thickEnd = end
        self.itemRgb = 0
        self.blockCount = 0
        self.exon_site = {}
        self.gene_biotype = gene_biotype
        # self.blockSizes = {}
        # self.blockStarts = {}

    def add_exon_block(self, s, e):
        # self.blockCount += 1
        self.exon_site["{0}-{1}".format(s, e)] = 0

    def tostr(self, biotype2bed13):
        blockSizes = None
        blockStarts = None
        blockCount = len(self.exon_site)

        for i, v in self.exon_site.items():
            tmp = i.split("-")
            if blockSizes is None:
                blockSizes = str(int(tmp[1]) - int(tmp[0]))
                blockStarts = tmp[0]
            else:
                blockSizes = blockSizes + "," + "{0}".format(int(tmp[1]) - int(tmp[0]))
                blockStarts = blockStarts + "," + tmp[0]

        obj_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(
            self.chrom,
            self.start,
            self.end,
            self.name,
            self.score,
            self.strand,
            self.thickStart,
            self.thickEnd,
            self.itemRgb,
            blockCount,
            blockSizes,
            blockStarts,
        )
        if biotype2bed13:
            obj_str += "\t{0}".format(self.gene_biotype)
        return obj_str


def _bedcol9_2dic(bedcol9, split=";"):
    tmp = bedcol9.split(split)
    col9dic = {}
    for i in tmp:
        i = i.strip()
        ii = i.split(" ")
        if len(ii) < 2:
            continue
        val = ii[1].strip('"')
        col9dic[ii[0]] = val
    return col9dic


def gtf2bed(
    gtf, outfile, sep, gene_name_tag, gene_id_tag, gene_biotype_tag, biotype2bed13
):
    gene = None
    gtf_hand = open(gtf, "r")
    gtf2bed = open(outfile, "w")

    for line in gtf_hand:
        bedtab = line.rstrip().split("\t")
        if len(bedtab) < 9:
            continue

        # if re.match('chr', bedtab[0]) == None:
        #    continue

        # if re.match('chrM', bedtab[0]) != None:
        #    continue

        # print(bedtab)
        # print(gene)

        if bedtab[2] == "gene":
            if gene is not None:
                gtf2bed.write(gene.tostr(biotype2bed13) + "\n")
            Gene_name = "xx"
            Gene_biotype = "yy"

            col9dic = _bedcol9_2dic(bedtab[8], split=sep)
            if gene_name_tag in col9dic.keys():
                Gene_name = col9dic[gene_name_tag]
            elif gene_id_tag in col9dic.keys():
                Gene_name = col9dic[gene_id_tag]
            else:
                print(bedtab[8], " have no {0}, xx instead".format(gene_name_tag))

            if gene_biotype_tag in col9dic.keys():
                Gene_biotype = col9dic[gene_biotype_tag]
            else:
                print(bedtab[8], " have no {0}, yy instead".format(gene_biotype_tag))

            gene = _genebed12(
                bedtab[0], bedtab[3], bedtab[4], Gene_name, bedtab[6], Gene_biotype
            )

        elif bedtab[2] == "exon":
            if gene is not None:
                gene.add_exon_block(bedtab[3], bedtab[4])
        else:
            pass

    # 最后一行没有输出

    gtf_hand.close()
    gtf2bed.close()
    print("gtf2bed finished")
