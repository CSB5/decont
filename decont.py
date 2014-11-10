#!/usr/bin/env python
"""Performs a mapping of given SR/PE reads (gzip supported) with
BWA-MEM against given source of contamination and produces an
(unsorted) BAM file with contaminated reads (one mate mapping suffices
to make pair count as contamination) and new, gzipped fastq file(s)
with uncontaminated reads.

Needs samtools and BWA(-MEM) installed.
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"


#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
from collections import deque
from collections import defaultdict
from string import maketrans
import gzip

#--- third-party imports
#
#/


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

"""
FIXME give flags a name

from samtools' 0.1.19 bam.h:

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
"""


def read_base_name(r):
    """return base name for read, i.e. without pairing information"""

    if r[-1] in "0123456789" and r[-2] in "#/":
        return r[:-2]
    else:
        return r


def complement(strand):
    """return DNA complement
    
    from http://stackoverflow.com/questions/1738633/more-pythonic-way-to-find-a-complementary-dna-strand
    """
    return strand.translate(maketrans('TAGCtagc', 'ATCGATCG'))


def sam_to_fastq(sam_line, fastq_fh, check_uniq_occurance=10000):
    """convert sam line to fastq entry

    Will make an attempt to check that no reads with identical names
    were processed (keeping check_uniq_occurance entries). The larger
    the number the more things will slow down. This makes it in fact
    the most time consuming step, even though it's in theory unnecessary
    as long as we pass 0x900==0 reads in here.
    """

    # local static fake
    # http://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
    if 'seen' not in sam_to_fastq.__dict__:
        sam_to_fastq.seen = deque(maxlen=check_uniq_occurance)

    line_split = sam_line.split('\t')
    (name, flag) = (line_split[0:2])
    name = read_base_name(name)# necessary?
    flag = int(flag)
    (seq, qual) = line_split[9:11]

    if flag & 0x10:
        seq = reversed(complement(seq))
        qual = reversed(qual)

    # From SAM spec: "If 0x1 is unset, no assumptions can be made
    # about 0x2, 0x8, 0x20, 0x40 and 0x80.". So check for pair before
    # checking read number.
    index = 1
    if (flag & 0x1) and (flag & 0x80):
        index += 1
    name = "%s/%d" % (name, index)

    if check_uniq_occurance:
        assert name not in sam_to_fastq.seen
        sam_to_fastq.seen.append(name)

    fastq_fh.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))


def main(fastq_in, ref, fastq_fh, bam_fh, num_threads=2, bwa='bwa'):
    """main function

    fastq_in and fastq_fh should be lists (single or paired-end)
    """

    # FIXME bufsize needs testing & optimization
    bufsize = 4096

    assert len(fastq_in) == len(fastq_fh)
    for f in fastq_in + [ref]:
        assert os.path.exists(f)

    # FIXME check support for bwa mem
    
    # can't use pysam with subprocess
    #p = subprocess.Popen(bwamem_cmd, stdout=subprocess.PIPE)
    #s = pysam.Samfile(p.stdout, "rb")
    # throws 'TypeError: Argument must be string or unicode'
    # see also https://www.biostars.org/p/15298/#105456
    # could use stringio?
    cmd = ["samtools", "view", "-bS", "-"]
    samtools_p = subprocess.Popen(cmd,
        stdin=subprocess.PIPE, stdout=bam_fh, bufsize=bufsize)

    # could redirect stderr to file and cat on problem but leave it to
    # user
    cmd = [bwa, 'mem', '-t', "%d" % num_threads, ref]
    cmd.extend(fastq_in)
    bwa_p = subprocess.Popen(cmd,
        stdout=subprocess.PIPE, bufsize=bufsize)

    if len(fastq_in) == 2:
        last_mate_name = None# for mate pair consitency check
    counts = defaultdict(int)# for reporting

    for line in bwa_p.stdout:
    #for (line_no, line) in enumerate(bwa_p.stdout):
        #if line_no>1000:
        #    LOG.critical("DEBUG break")
        #    break

        # catch SAM header
        if line.startswith('@'):
            samtools_p.stdin.write(line)
            continue

        counts['received from BWA'] += 1
        (name, flag) = line.split('\t')[:2]
        flag = int(flag)
        name = read_base_name(name)# necessary?

        # from the SAM spec:
        # http://samtools.github.io/hts-specs/SAMv1.pdf:
        #
        # - Bit 0x4 is the only reliable place to tell whether the
        # read is unmapped. If 0x4 is set, no assumptions can be made
        # about RNAME, POS, CIGAR, MAPQ,
        #
        # - If 0x1 is unset, no assumptions can be made about 0x2,
        # 0x8, 0x20, 0x40 and 0x80.
        #
        # - For each read/contig in a SAM file, it is required that
        # one and only one line associated with the read satisfies
        # 'FLAG & 0x900 == 0'. This line is called the primary line of
        # the read.
        #
        # - Bit 0x100 marks the alignment not to be used in certain
        # analyses when the tools in use are aware of this bit. It is
        # typically used to flag alternative mappings when multiple
        # mappings are presented in a SAM.
        #
        # - Bit 0x800 indicates that the corresponding alignment line
        # is part of a chimeric alignment. A line flagged with 0x800 is
        # called as a supplementary line

        if flag & 0x900:
            LOG.debug("Skipping secondary or"
                      " supplementary alignment: %s" % line)
            continue

        is_paired = flag & 0x1
        if len(args.fastq_in) == 2:
            assert is_paired, (
                "Got two fastq files (paired end) but bwa mem output"
                 " reported single end mapping")
        else:
            assert not is_paired, (
                "Got one fastq files (single end) but bwa mem output"
                 " reported paired end mapping")

        if is_paired:
            if last_mate_name:
                assert name == last_mate_name
                last_mate_name = None
            else:
                last_mate_name = name

            both_unmapped = (flag & 0x4) and (flag & 0x8)
            if not both_unmapped:
                counts['written to %s' % bam_fh.name] += 1
                samtools_p.stdin.write(line)
            else:
                if flag & 0x40:# 1st in pair
                    counts['written to %s' % fastq_fh[0].name] += 1
                    sam_to_fastq(line, fastq_fh[0])
                elif flag & 0x80:# 2nd in pair
                    counts['written to %s' % fastq_fh[1].name] += 1
                    sam_to_fastq(line, fastq_fh[1])
                else:
                    raise ValueError()
        else:
            unmapped = flag & 0x4
            if unmapped:
                counts['written to %s' % fastq_fh[0].name] += 1
                sam_to_fastq(line, fastq_fh[0])
            else:
                counts['written to %s' % bam_fh.name] += 1
                samtools_p.stdin.write(line)

    samtools_p.stdin.close()
    if samtools_p.wait() != 0:
        LOG.critical("Unhandled samtools error. Note, this can happen"
                     " when BAM is empty i.e. with no contamination."
                     " samtools might then issue the following complaint:"
                     " \"reference 'XYZ' is recognized as '*'\""
                     " followed by \"truncated file\"")

    for (k, v) in counts.items():
        LOG.info("Reads %s: %d" % (k, v))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # FIXME support bwa-mem special args
    # FIXME allow reuse of existing bam (needs to be sorted by name!)
    mandatory = parser.add_argument_group('mandatory arguments')
    mandatory.add_argument("-i", "--fq",
                        dest='fastq_in',
                        required=True,
                        nargs="*",
                        help="FastQ input file/s")
    mandatory.add_argument("-o", "--outpref",
                        dest='outpref',
                        required=True,
                        help="Filename prefix for output files")
    mandatory.add_argument("-r", "--ref",
                        required=True,
                        dest='ref',
                        help="Reference fasta file of source of contamination"
                        " (needs to be bwa indexed already)")
    
    default = 8
    parser.add_argument("-t", "--threads",
                        dest='num_threads',
                        default=default,
                        type=int,
                        help="Number of threads to use for mapping"
                        " (default = %d)" % default)
    parser.add_argument("-b", "--bwa",
                        dest='bwa',
                        help="Path to BWA supporting the mem command (will use BWA found in PATH if not set)")
    args = parser.parse_args()

    fastq_out = ["%s_%d.fastq.gz" % (args.outpref, i+1)
                 for i in range(len(args.fastq_in))]
    bam_out = "%s.bam" % (args.outpref)

    infiles = args.fastq_in + [args.ref]
    for f in infiles:
        if not os.path.exists(f):
            LOG.fatal("Input file %s does not exist" % f)
            sys.exit(1)

    outfiles = fastq_out + [bam_out]
    for f in outfiles:
        if os.path.exists(f):
            LOG.fatal("Cowardly refusing to overwrite"
                      " already existing file %s" % f)
            sys.exit(1)

    fastq_fh = [gzip.open(f, 'w') for f in fastq_out]

    bam_fh = open(bam_out, 'wb')

    if not os.path.exists(args.ref + ".bwt"):
        LOG.warn("Doesn't look like reference was indexed."
                 " Did you forget to run bwa index %s ?" % args.ref)

    main(args.fastq_in, args.ref, fastq_fh, bam_fh,
         num_threads=args.num_threads, bwa=args.bwa)

    for f in fastq_fh + [bam_fh]:
        f.close()

    LOG.info("Successful program exit")
