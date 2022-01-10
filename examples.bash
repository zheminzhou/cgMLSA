python modules/01_alignAlleles.py -p examples/examples.profile.gz -n 10 --ublast -d examples/genes/ examples.examples.alleles.fasta.gz
python modules/02B_alleleAlign2GenomeDistance.py -o examples/examples.dist -p examples/examples.profile.gz -d examples/genes
python modules/03_getAlleleTrees.py -d examples/genes --raxml_ng
python modules/04_expendAlleleTree.py -d examples/genes -p examples/examples.profile.gz -o examples/examples.geneTrees
python modules/05_summariseGeneTrees.py -i examples/examples.geneTrees -o examples/examples.guide_tree
python modules/06_splitAstridTreeByQuartetSupport.py -i examples/examples.guide_tree -o examples/examples.subset -m 10 -x 20
python modules/07_generateSubtreesForEachGroup.py -g examples/examples.subset -i examples/examples.geneTrees -p examples/examples
python modules/08_runAstralForEachSubtrees.py -d examples
python modules/09_summarizeAstralSubtrees2SuperTree.py -d examples -o examples/examples.supertree
python modules/10_reviseBranchLength.py -t examples/examples.supertree -g examples/examples.geneTrees -d examples/examples.dist -o examples/examples.species_tree
