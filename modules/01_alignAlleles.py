import os, click, numpy as np, subprocess, tempfile, pandas as pd, re
from uberBlast import uberBlast
from configure import readFasta, uopen, rc, externals
from multiprocessing import Pool

mafft = externals['mafft']
def alignAGene(data) :
    fin, fout = data
    subprocess.Popen('{0} --thread 1 --auto {1} | gzip > {2}'.format(mafft, fin, fout), stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()
    return  fout


def ublast_get_seq(seq, align) :
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

def alignAGene_ublast(data) :
    infile, outfile = data

    seqs = readFasta(infile)
    with tempfile.TemporaryDirectory(dir='.', prefix='ub_') as tmpdir :
        with open(os.path.join(tmpdir, 'ref'), 'wt') as fout :
            for n, s in seqs.items() :
               fout.write('>{0}\n{1}\n'.format(n, s))
               break
        aligns = uberBlast(
            '-r {0} -q {1} -f -m --blastn --diamond --min_id 0.4 --min_cov 30 --min_ratio 0.2 --merge_gap 600 --merge_diff 1.5 -t 1 -s 1 -e 0,3'.format( \
                os.path.join(tmpdir, 'ref'), infile,
            ).split())

    out_seqs = {}
    for align in aligns:
        n = align[0]
        seq = ublast_get_seq(seqs[n], align)
        if n not in out_seqs:
            out_seqs[n] = seq
        else:
            s2 = ''.join([m if m != '-' else n for m, n in zip(out_seqs[n], seq)])
            out_seqs[n] = s2

    with uopen(outfile, 'w') as fout:
        for n, seq in out_seqs.items():
            fout.write('>{0}\n{1}\n'.format(n, seq))
    return outfile


def get_allele_names(profile_file) :
    profile = pd.read_csv(profile_file, sep='\t', header=0, dtype=str)
    genes = profile.columns
    idx = np.array([ not (i == 0 or g.startswith('#')) for i, g in enumerate(genes) ])
    genes = np.array([g.split(':')[-1] for g in genes[idx]], dtype=str)
    alleles = profile.values[:, idx].T.astype(int)

    allele_seqs = {}
    for g, ids in zip(genes, alleles) :
        allele_seqs[g] = {}
        for i in np.unique(ids[ids>0]) :
            allele_seqs[g]['{0}_{1}'.format(g, i)] = []
    return allele_seqs


@click.command()
@click.option('-p', '--profile', help='profile file (can in gzip format)')
@click.option('-d', '--outdir', help='folder storing output alignment')
@click.option('-n', '--n_proc', help='n_proc', default=10, type=int)
@click.option('--ublast', help='flag to use uBlast (EToKi) instead of mafft for alignment', default=False, is_flag=True)
@click.argument('alleles', nargs=-1)
def alignAlleles(profile, alleles, outdir, ublast, n_proc) :
    global pool
    pool = Pool(n_proc)
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    allele_seqs = get_allele_names(profile)

    for fname in alleles :
        with uopen(fname, 'rt') as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0].split(':')[-1]
                    gene, idx = name.rsplit('_', 1)
                elif name in allele_seqs.get(gene, {}) :
                    allele_seqs[gene][name].extend(line.split())
    allele_files = []
    for gene, seqs in allele_seqs.items() :
        allele_file = os.path.join(outdir, '{0}.alleles.fasta'.format(gene))
        aln_file = os.path.join(outdir, '{0}.aln.fasta.gz'.format(gene))
        allele_files.append([allele_file, aln_file])
        with open(allele_file, 'wt') as fout :
            for n, s in seqs.items():
                assert len(s) > 0, 'allele {0} does not have a sequence.'.format(n)
                fout.write('>{0}\n{1}\n'.format(n, ''.join(s)))
    outputs = []
    if ublast :
        for aln_seqs in pool.imap_unordered(alignAGene_ublast, allele_files) :
            outputs.append(aln_seqs)
    else :
        for aln_seqs in pool.imap_unordered(alignAGene, allele_files):
            outputs.append(aln_seqs)
    return outputs

pool = None
if __name__ == '__main__' :
    alignAlleles()
