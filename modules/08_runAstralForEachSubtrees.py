import subprocess

import click, os, glob
try :
    from .configure import externals
except :
    from configure import externals

astral = externals['astral']


@click.command()
@click.option('-d', '--folder', help='folder storing gene trees', default='.')
@click.option('-i', '--input_pattern', help='Default: *.__G_*__', default='*.__G_*__')
@click.option('-r', '--repeats', help='number of sampling per tree', default=2)
def main(folder, input_pattern, repeats):
    subtrees = glob.glob(os.path.join(folder, input_pattern))
    for subtree in subtrees :
        subprocess.Popen('{0} -i {1} -o {1}.astral -w {2} -t 3 -C -T 10'.format(astral, subtree, repeats).split()).communicate()



if __name__ == '__main__':
    main()
