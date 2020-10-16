import click, subprocess
from multiprocessing import Pool
try :
    from .configure import externals, uopen
except :
    from configure import externals, uopen

fasttree = externals['fasttree']

def run_tree(alignment) :
    subprocess.Popen('gzip -cd {1}|{0} -gtr -pseudo -nt -quiet > {2}'.format(fasttree, alignment, alignment+'.'+pf), stdout=subprocess.PIPE, shell=True).communicate()
    

@click.command()
@click.option('-n', '--n_proc', help='number of processes. Default: 10', default=10, type=int)
@click.option('-p', '--postfix', help='postfix added to the outputfile. default: fasttree', default='fasttree')
@click.argument('alignments', nargs=-1)
def alignAlleles(n_proc, postfix, alignments) :
    global pf
    pf = postfix
    pool = Pool(n_proc)
    list(pool.map(run_tree, alignments))

pf = None
if __name__ == '__main__' :
    alignAlleles()
