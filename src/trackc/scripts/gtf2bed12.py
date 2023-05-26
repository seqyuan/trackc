import argparse

class _genebed12(object):
    """bed12 object"""
    
    def __init__(self,chrom,start,end,name,strand, gene_biotype):
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
        #self.blockSizes = {}
        #self.blockStarts = {}

    def add_exon_block(self,s,e):        
        #self.blockCount += 1
        self.exon_site["{0}-{1}".format(s,e)] = 0
        
            
    def tostr(self, biotype2bed13):
        blockSizes = None
        blockStarts = None
        blockCount = len(self.exon_site)
        
        for i,v in self.exon_site.items():
            tmp = i.split("-")
            if blockSizes == None:
                blockSizes = str(int(tmp[1]) - int(tmp[0]))
                blockStarts = tmp[0]
            else:
                blockSizes = blockSizes + "," + "{0}".format(int(tmp[1]) - int(tmp[0]))
                blockStarts = blockStarts + "," + tmp[0]
        
      
        obj_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}".format(self.chrom, self.start, self.end,\
                                                                                        self.name, self.score, self.strand,\
                                                                                        self.thickStart, self.thickEnd, self.itemRgb,\
                                                                                blockCount, blockSizes, blockStarts, self.gene_biotype)
        if biotype2bed13 == False:
             obj_str = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(self.chrom, self.start, self.end,\
                                                                                        self.name, self.score, self.strand,\
                                                                                        self.thickStart, self.thickEnd, self.itemRgb,\
                                                                                blockCount, blockSizes, blockStarts)
        return(obj_str)
            

def _bedcol9_2dic(bedcol9, split=";"):
    tmp = bedcol9.split(split)
    col9dic = {}
    for i in tmp:
        i = i.strip()
        ii = i.split(" ")
        if len(ii) < 2:
            continue
        val = ii[1].strip("\"")
        col9dic[ii[0]] = val
    return col9dic



gene = None
#gene_id="gene_id"
gene_name_tag="gene_name"
gene_id_tag="gene_id"
gene_biotype_tag="gene_biotype"
split = ";"
biotype2bed13 = True


def _parse_arguments(args=None):
    """
    get command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
        description='Converts the GTF to bed12 or bed13')

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--gtf', '-g',
                                help='input file. GTF format gene annotation.',
                                required=True)

    parserRequired.add_argument('--outFileName', '-o',
                                help='File name of gene annotation bed12 format.',
                                required=True)

    parserOpt = parser.add_argument_group('Optional arguments')
    parserOpt.add_argument('--splitby', '-s',
                           help='tags of GTF split by which str',
                           default=";")
    parserOpt.add_argument('--gene_name_tag', '-nt',
                           help='tags of GTF gene_name',
                           default="gene_name")
    parserOpt.add_argument('--gene_id_tag', '-it',
                           help='tags of GTF gene_id',
                           default="gene_id")
    parserOpt.add_argument('--gene_biotype_tag', '-bt',
                           help='tags of GTF gene_biotype',
                           default="gene_biotype")
    parserOpt.add_argument('--biotype2bed13', '-bed13',
                           help='if this arg is set, the out file columns 13 will add gene\'s biotype',
                           action='store_true')
    parserOpt.add_argument("-h", "--help", action="help", help="show this help message and exit")

    return parser


def _main(args=None):

    args = _parse_arguments().parse_args(args)
    gene_name_tag = args.gene_name_tag
    gene_id_tag = args.gene_id_tag
    gene_biotype_tag = args.gene_biotype_tag
    split = args.splitby
    biotype2bed13 = args.outFileName

    gene = None
    gtf_hand = open(args.gtf, "r")
    gtf2bed = open(args.outFileName,'w')

    for line in gtf_hand:
        bedtab = line.rstrip().split('\t')
        if len(bedtab) < 9:
            continue

        #if re.match('chr', bedtab[0]) == None:
        #    continue
            
        #if re.match('chrM', bedtab[0]) != None:
        #    continue
            
        #print(bedtab)
        #print(gene)
            
        if bedtab[2] == "gene":
            if gene != None:
                gtf2bed.write(gene.tostr(biotype2bed13) + "\n")
            Gene_name = 'xx'
            Gene_biotype = "yy"
            
            
            col9dic = _bedcol9_2dic(bedtab[8], split=split)
            if gene_name_tag in col9dic.keys():
                Gene_name = col9dic[gene_name_tag]
            elif gene_id_tag in col9dic.keys():
                Gene_name = col9dic[gene_id_tag]
            else:
                print(bedtab[8], ' have no {0}, xx instead'.format(gene_name_tag))
            
            if gene_biotype_tag in col9dic.keys():
                Gene_biotype = col9dic[gene_biotype_tag]
            else:
                print(bedtab[8], ' have no {0}, yy instead'.format(gene_biotype_tag))
                
            gene = _genebed12(bedtab[0], bedtab[3], bedtab[4], Gene_name, bedtab[6], Gene_biotype)
        
        elif bedtab[2] == "exon":
            if gene != None:
                gene.add_exon_block(bedtab[3],bedtab[4])
        else:
            pass
            
    #最后一行没有输出
            
    gtf_hand.close()
    gtf2bed.close()
    print ('gtf2bed finished')
