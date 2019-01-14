PGHM, Version 1.0

Installation and Usage


INSTALLATION

python setup.py install


DEPENDENCIES

matplotlib
NumPy




USAGE

PGHM is intended to be easy to use. For typical usage for comparing two samples, PGHM requires 4 input parameters:

1) Sequences file from two samples, each sequence classified into a phylogenetic group (see test.fasta file for format)
2) Accession Number of sample 1
3) Accession Number of sample 2
3) Output phylogenetic heatmap figure file 

The sequences from the two samples are classified into four phylogenetic groups and combined together in one file.
The file is sorted.

Example using the test data,

./PGHM.py -i test.fasta -1 SRR10000000 -2 SRR20000000 -o figure

## Figure
![figure](https://user-images.githubusercontent.com/44556007/51094474-553d2e80-177b-11e9-94a6-2db57f0cc72b.png)
