# nextstrain-ToBRFV

In October 2019, the Dutch national reference centre for plant health received a tomato sample suspected to be infected with ToBRFV. Presence of the virus was confirmed with analysis of plant rRNA depleted Illumina RNAseq data, using our virus detection pipeline. A trace-back survey was initiated to identify possible linkages between ToBRFV genotypes, companies and epidemiological traits. Nextstrain was used to visualize these potential connections. 

Nextstrain implementation

Nextstrain is a bioinformatics pipeline that uses two tools, Augur and auspice. Augur (github.com/Nextstrain/augur) is a bioinformatics toolkit for phylogenetic analysis that uses a series of modules that produces output that subsequently can be visualized by auspice (github.com/Nextstrain/auspice). 

All files necessary to run Augur are present in the data and config directory. In the data directory a multi-fasta file contains all included ToBRFV genomes and all metadata that belong to these genomes can be found in the .tsv file. The config directory contains the lat-long file, the auspice config file, and the annotated reference genome file, which is necessary to determine the amino acid changes. Additional features had to be added to the NCBI Genbank file (KT383474_NCBI-edit.gb) to enable the use of this sequence as reference, i.e. gene names and translation tables for the four CDS annotations.

We included the two json files that are necessary for visualization by auspice. However, these files can also be generated using the Snakefile included, which runs each step in Augur. Before running the Snakefile, empty (or rename the json files in) the auspice directory and add an empty results directory.
