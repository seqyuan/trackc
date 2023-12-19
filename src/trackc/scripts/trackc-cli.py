import click

@click.command()
@click.option('--regions', '-r', required=True, 
              help="genome regions, eg.: 18:47950000-48280000 18:75280000-74850000")
@click.option('--outfile', '-o', default='trackc.pdf',
              help='output filename')
def hello(regions, outfile):
    print(1111, regions)
    print(1111, outfile)
    

if __name__ == '__main__':
    hello()




