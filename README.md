decont
======

Decontaminate FastQ files by mapping with BWA-MEM against a given
source.


Performs a mapping of given SR/PE reads (gzip supported) with BWA-MEM
against given source of contamination and produces an (unsorted) BAM
file with contaminated reads (one mate mapping suffices to make pair
count as contamination) and new, gzipped fastq file(s) with
uncontaminated reads.

Needs samtools and BWA(-MEM) installed.

