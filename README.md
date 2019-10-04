decont
======


# NOTE: this repo is outdated. Please use the fully reimplemented and faster version called [decontanimate](https://bitbucket.org/andreas-wilm/decontanimate/src/master/)

Decontaminate FastQ files by mapping with BWA-MEM against a given
source.


Performs a mapping of given SR/PE reads (gzip supported) with BWA-MEM
against given source of contamination and produces an (unsorted) BAM
file with contaminated reads (one mate mapping suffices to make pair
count as contamination) and new, gzipped fastq file(s) with
uncontaminated reads.

Needs samtools and BWA(-MEM) installed.

If minimum coverage is not met, the read is treated as unaligned.

# TODO


## Add testing data

- SE aligned
- SE unaligned
- SE secondary
- SE aligned below mincov
- PE both aligned
- PE both unaligned
- PE one aligned
- PE both aligned secondary



## Remove most python bits for speedup

- Run `samtools fixmate` (paranoia) followed by
- `samtools fastq -F`, removing any unmapped read or read with unmapped
pair (option) and secondary alignments
