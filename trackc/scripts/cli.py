#!/Users/yuanzan/anaconda3/bin/python3
# -*- coding: utf-8 -*-

import warnings

import click

warnings.filterwarnings("ignore")

import trackc as tc
from trackc.scripts.gtf2bed12 import gtf2bed
from trackc.scripts.trackc_cli import cli as track

# from trackc_cli import track
# from trackc.scripts.tadscore import tadScore


@click.group()
def main():
    pass


@main.command(
    name="cli",
    short_help="plot trackc by use yaml format configfile, yaml track type parameters is the same as API",
)
@click.argument("config", metavar="<trackc-conf.yml>")
@click.option(
    "--regions",
    "-r",
    default=None,
    help="genome regions, \
              eg.:'18:47950000-48280000 18:75280000-74850000'",
)
@click.option(
    "--outfile", "-o", default="trackc.pdf", help="output filename. default=trackc.pdf"
)
@click.option(
    "--basefigsize",
    "-s",
    default="6,1",
    help='base figsize: width,height. default="6,1"\
              The height option in the config.yml is relative to the base figsize height',
)
def cli(config, regions, outfile, basefigsize):
    track(config, regions, outfile, basefigsize)


@main.command(
    name="gtf2bed", short_help="Converts the gene GTF to bed12 or bed13 format"
)
@click.argument("gtf", metavar="<genome.gtf>")
# , help="input file. GTF format gene annotation."
@click.option(
    "--outfile",
    "-o",
    default="gene.bed",
    help="File name of gene annotation bed12 format.",
)
@click.option("--sep", "-s", default=";", help="tags of GTF sepration by which str")
@click.option(
    "--gene_name_tag", "-nt", default="gene_name", help="tags of GTF gene_name"
)
@click.option("--gene_id_tag", "-it", default="gene_id", help="tags of GTF gene_id")
@click.option(
    "--gene_biotype_tag",
    "-bt",
    default="gene_biotype",
    help="tags of GTF gene_biotype (default: gene_biotype, fallback: gene_type)",
)
@click.option(
    "--biotype2bed13/--no-biotype2bed13",
    "-bed13",
    is_flag=True,
    show_default=True,
    default=False,
    help="if this arg is set, the out file columns 13 will add gene's biotype",
)
def mode2(
    gtf, outfile, sep, gene_name_tag, gene_id_tag, gene_biotype_tag, biotype2bed13
):
    gtf2bed(
        gtf, outfile, sep, gene_name_tag, gene_id_tag, gene_biotype_tag, biotype2bed13
    )


# @main.command()
# @click.argument('mcool', metavar='name.mcool::/resolutions/10000')

# def tadscore(mcool):
#    tadScore()


if __name__ == "__main__":
    main()
