# cgMLSA
phylogenetic analytic toolsets based on cgMLST

# INSTALLATION:

cgMLSA was developed and tested in Python 3.8. cgMLSA depends on several Python libraries: 
~~~~~~~~~~
ete3
numba
numpy
pandas
click
~~~~~~~~~~

All libraries can be installed using pip: 

~~~~~~~~~~
pip install ete3 numba numpy pandas click
~~~~~~~~~~
EToKi also calls the following 3rd party programs for different pipelines:

~~~~~~~~~~
Astral (v5.15.1)
ASTRID (v1.4)
blastn
diamond
erable
FastTreeMP
mafft
makeblastdb
rapidnj
~~~~~~~~~~

All 3rd party programs need to be installed in the cgMLST/dependencies folder. 



# Quick Start (with examples)

### Build multi-sequence alignments for every genes in the cgMLST scheme
~~~~~~~~~~~
python modules/01_alignAlleles.py -p examples/examples.profile.gz -n 10 --ublast -d examples/genes/ examples.examples.alleles.fasta.gz
~~~~~~~~~~~
### Calculate nucleotide distances of genomes
~~~~~~~~~~~
python modules/02B_alleleAlign2GenomeDistance.py -o examples/examples.dist -p examples/examples.profile.gz -d examples/genes
~~~~~~~~~~~
### Calculate allelic trees
~~~~~~~~~~~
python modules/03_getAlleleTrees.py -d examples/genes --raxml_ng
~~~~~~~~~~~
### Add genomic labels into the gene trees
~~~~~~~~~~~
python modules/04_expendAlleleTree.py -d examples/genes -p examples/examples.profile.gz -o examples/examples.geneTrees
~~~~~~~~~~~
### Calculate a guide_tree using ASTRID, and calculate posterior probabilities (PP) of branches using ASTRAL
~~~~~~~~~~~
python modules/05_summariseGeneTrees.py -i examples/examples.geneTrees -o examples/examples.guide_tree
~~~~~~~~~~~
### Split the guide tree based on PP values
~~~~~~~~~~~
python modules/06_splitAstridTreeByQuartetSupport.py -i examples/examples.guide_tree -o examples/examples.subset -m 10 -x 20
~~~~~~~~~~~
### Generate disjoint subsets of gene trees
~~~~~~~~~~~
python modules/07_generateSubtreesForEachGroup.py -g examples/examples.subset -i examples/examples.geneTrees -p examples/examples
~~~~~~~~~~~
### Calculate ASTRAL trees for each disjoint subsets
~~~~~~~~~~~
python modules/08_runAstralForEachSubtrees.py -d examples
~~~~~~~~~~~
### Generate the supertree by concatenating disjoint subtrees together
~~~~~~~~~~~
python modules/09_summarizeAstralSubtrees2SuperTree.py -d examples -o examples/examples.supertree
~~~~~~~~~~~
### Update the branch lengths of the supertree based on genomic distances
~~~~~~~~~~~
python modules/10_reviseBranchLength.py -t examples/examples.supertree -g examples/examples.geneTrees -d examples/examples.dist -o examples/examples.species_tree
~~~~~~~~~~~
