import click

@click.command()
@click.argument('--regions', '-r', required=True, 
              help="genome regions, eg.: 18:47950000-48280000 18:75280000-74850000")
@click.argument('--outfile', '-o', default='trackc.pdf',
              help='output filename')

@click.option('--baseFigsize', '-s', default='8,1',
              help='base figsize: width,height . \
              The height option in the config.yml is relative to the base figsize height')


def hello(regions, outfile):
    print(1111, regions)
    print(1111, outfile)
    

if __name__ == '__main__':
    hello()




