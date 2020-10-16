import os, sys, click, numpy as np, pandas as pd, re
from uberBlast import uberBlast
from configure import readFasta, uopen, rc
from _collections import OrderedDict

def get_seq(seq, align) :
    seq2 = seq[align[6]-1:align[7]]
    if align[8] < align[9] :
        s, e = align[8:10]
    else :
        s, e = align[9], align[8]
        seq2 = rc(seq2)
    out_seq = []
    idx = 0
    for r, t in re.findall(r'(\d+)([MDI])', align[14]) :
        r = int(r)
        if t == 'M' :
            out_seq.append(seq2[idx:idx+r])
            idx +=r
        elif t == 'D' :
            out_seq.append('-'*r)
        else :
            idx += r
    return ''.join(['-'*(s-1)] + out_seq + ['-'*(align[13]-e)])

def processGenes(data, outdir, n_proc) :
    for gene, seqs in data.items() :
        with uopen('ref', label='x') as reffile, uopen('qry', label='x') as qryfile:
            for n, s in seqs.items():
                seqs[n] = ''.join(s)
                qryfile.write('>{0}\n{1}\n'.format(n, seqs[n]))
            qryfile.close()
            r = list(seqs.items())[0]
            reffile.write('>{0}\n{1}\n'.format(r[0], r[1]))
            reffile.close()

            aligns = uberBlast(
                '-r {0} -q {1} -f -m --blastn --diamond --min_id 0.4 --min_cov 30 --min_ratio 0.2 --merge_gap 600 --merge_diff 1.5 -t {2} -s 1 -e 0,3'.format( \
                    reffile.name, qryfile.name, n_proc
                ).split())

        out_seqs = {}
        for align in aligns:
            n = align[0]
            seq = get_seq(seqs[n], align)
            if n not in out_seqs:
                out_seqs[n] = seq
            else:
                s2 = ''.join([m if m != '-' else n for m, n in zip(out_seqs[n], seq)])
                out_seqs[n] = s2

        fout = os.path.join(outdir, '{0}.aln.fasta.gz'.format(gene))
        with uopen(fout, 'w') as fout:
            for n, seq in out_seqs.items():
                fout.write('>{0}\n{1}\n'.format(n, seq))


@click.command()
@click.option('-i', '--input', help='allelic sequences (can in gzip format)')
@click.option('-o', '--outdir', help='folder storing output alignment')
@click.option('-p', '--n_proc', help='n_proc', default=10, type=int)
def alignGeneSeqs(input, outdir, n_proc) :
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with uopen(input, 'r') as fin:
        data = {}
        for line in fin:
            if line.startswith('>'):
                name = line[1:].strip().split(' ', 1)[0]
                gene, allele = name.rsplit('_', 1)
                if gene not in data:
                    processGenes(data, outdir, n_proc)
                    data = {gene: OrderedDict() }
                data[gene][name] = []
            else :
                data[gene][name].extend(line.strip().split())
    processGenes(data, outdir, n_proc)

if __name__ == '__main__' :
    alignGeneSeqs()