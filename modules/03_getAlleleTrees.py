import shutil
import click, subprocess, tempfile, gzip, os, glob
from multiprocessing import Pool
try :
    from .configure import externals, uopen
except :
    from configure import externals, uopen

fasttree = externals['fasttree']
raxml = externals['raxml_ng']

def run_tree(data) :
    fn, raxml_ng = data
    with tempfile.TemporaryDirectory(dir='.', prefix='at_') as tmpdir :
        f = os.path.join(tmpdir, 'aln.fas')
        with gzip.open(fn, 'rt') as fin, open(f, 'wt') as fout :
            fout.write(fin.read())
        if raxml_ng :
            cmd = '{0} --thread 1 --redo --force --msa {1} --precision 8 --model GTR+G --blmin 1e-7 --site-repeats on --tree pars{{2}}'.format(raxml, f)
            subprocess.Popen(cmd.split(), stdout=subprocess.PIPE).communicate()
            if os.path.isfile(f+'.raxml.bestTree') :
                shutil.copy(f+'.raxml.bestTree', fn+'.'+pf)
        else :
            subprocess.Popen('{0} -gtr -pseudo -nt -quiet > {2}'.format(fasttree, f, fn+'.'+pf), stdout=subprocess.PIPE, shell=True).communicate()
    

@click.command()
@click.option('-n', '--n_proc', help='number of processes. Default: 10', default=10, type=int)
@click.option('-i', '--input_pattern', help='Default: *.aln.fasta.gz', default='*.aln.fasta.gz')
@click.option('-p', '--postfix', help='postfix added to the outputfile. default: tree', default='tree')
@click.option('-r', '--raxml_ng', help='flag to user raxml_ng instead of fasttree [default]', default=False, is_flag=True)
@click.option('-d', '--dir', help='folder storing output alignment')
def alignAlleles(n_proc, postfix, dir, input_pattern, raxml_ng) :
    alignments = glob.glob(os.path.join(dir, input_pattern))
    global pf
    pf = postfix
    pool = Pool(n_proc)
    list(pool.map(run_tree, [ [fn, raxml_ng] for fn in alignments]))

pf = None
if __name__ == '__main__' :
    alignAlleles()
