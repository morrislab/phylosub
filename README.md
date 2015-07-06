Deprecation notice
==================
PhyloSub has been replaced by
[PhyloWGS](http://www.genomebiology.com/2015/16/1/35), which performs the same
function as PhyloSub but with numerous improvements. Prime amongst these is the
ability to integrate copy number variation (CNV) data. PhyloWGS can, however,
run without any CNV information -- simply run with an empty `cnv_data.txt`, and
PhyloWGS will produce phylogenies based solely on single-nucleotide somatic
mutation (SSM) frequency, just as PhyloSub would.

This Python/C++ code is the accompanying software for the paper:

[Wei Jiao, Shankar Vembu, Amit G. Deshwar, Lincoln Stein and Quaid Morris, Inferring clonal evolution of tumors from single nucleotide somatic mutations, BMC Bioinformatics 15:35, 2014.](http://www.biomedcentral.com/1471-2105/15/35)

Please check the [README file in c++ folder](c++/README) for instructions to run the software.
