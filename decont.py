#!/usr/bin/env python
"""Performs a mapping of given SR/PE reads (gzip supported) with
BWA-MEM against given source of contamination and produces an
(unsorted) BAM file with contaminated reads (one mate mapping suffices
to make pair count as contamination) and new, gzipped fastq file(s)
with uncontaminated reads.

Needs samtools and BWA(-MEM) installed.
"""



#--- standard library imports
#
import sys
import logging
import os
import argparse
import subprocess
#from collections import deque
from collections import defaultdict
from string import maketrans
import gzip
from itertools import groupby
from collections import namedtuple
import tempfile


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014-2016 Genome Institute of Singapore"
__license__ = "WTFPL http://www.wtfpl.net/"



SamRead = namedtuple('SamRead',
                     ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar',
                      'rnext', 'pnext', 'tlen', 'seq', 'qual', 'rest'])

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



def read_add_index(r):
    """index is part of filename in fastq but not in BAM"""

    if r.qname == read_base_name(r.qname):
        # From SAM spec: "If 0x1 is unset, no assumptions can be made
        # about 0x2, 0x8, 0x20, 0x40 and 0x80.". So check for pair before
        # checking read number.
        index = 1
        if (r.flag & 0x1) and (r.flag & 0x80):
            index += 1
        r = r._replace(qname="%s/%d" % (r.qname, index))

    return r


def sam_to_read(line):
    """converts same line to namedtuple read
    """

    ls = line.rstrip().split('\t')
    l = ls[:11]
    l.append('\t'.join(ls[11:]))
    r = SamRead._make(l)
    r = r._replace(flag=int(r.flag))
    r = r._replace(pos=int(r.pos))
    r = r._replace(mapq=int(r.mapq))
    r = r._replace(pnext=int(r.pnext))
    r = r._replace(tlen=int(r.tlen))

    return r


def read_to_sam(read):
    """convert namedtuple read to sam line"""

    values = read._asdict().values()
    sam = '\t'.join([str(v) for v in values[:11]])
    if read.rest:
        sam = "%s\t%s" % (sam, read.rest)
    return sam + "\n"


def complement(strand):
    """return DNA complement

    from http://stackoverflow.com/questions/1738633/more-pythonic-way-to-find-a-complementary-dna-strand

    >>> complement('AcGtN')
    'TGCAN'
    """
    return strand.translate(maketrans('TAGCtagc', 'ATCGATCG'))


def read_to_fastq(read, fastq_fh):
    """convert sam line to fastq entry

    Will make an attempt to check that no reads with identical names
    were processed (keeping check_uniq_occurance entries). The larger
    the number the more things will slow down. This makes it in fact
    the most time consuming step, even though it's in theory unnecessary
    as long as we pass 0x900==0 reads in here.
    """

    # local static fake
    # http://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
    #if 'seen' not in sam_to_fastq.__dict__:
    #    sam_to_fastq.seen = deque(maxlen=check_uniq_occurance)

    # not part of SAM but part of fastq
    read = read_add_index(read)

    if read.flag & 0x10:
        read = read._replace(seq=complement(read.seq)[::-1])
        read = read._replace(qual=read.qual[::-1])

    #if check_uniq_occurance:
    #    assert name not in sam_to_fastq.seen
    #    sam_to_fastq.seen.append(name)

    fastq_fh.write('@%s\n%s\n+\n%s\n' % (read.qname, read.seq, read.qual))


def read_base_name(qname):
    """return base name for read, i.e. without pairing information

    >>> r = "M01853:160:000000000-ADF3H:1:1101:15677:1332/1"
    >>> read_base_name(r)
    'M01853:160:000000000-ADF3H:1:1101:15677:1332'

    >>> r = "M01853:160:000000000-ADF3H:1:1101:15677:1332#2"
    >>> read_base_name(r)
    'M01853:160:000000000-ADF3H:1:1101:15677:1332'

    >>> r = "M01853:160:000000000-ADF3H:1:1101:15677:1332"
    >>> read_base_name(r)
    'M01853:160:000000000-ADF3H:1:1101:15677:1332'
    """

    if qname[-1] in "0123456789" and qname[-2] in "#/":
        return qname[:-2]
    else:
        return qname


def bwa_mem_support(bwa='bwa'):
    """checks whether bwa binary supports bwa mem"""

    p = subprocess.Popen(bwa, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stderr:
        if 'mem' in line.split():
            return True
    return False


def cigar_items(cigar_str):
    """Taken from Brent Pedersen's https://github.com/brentp/cigar/blob/master/cigar.py
    """
    if cigar_str == "*":
        yield (0, None)
        raise StopIteration
    cig_iter = groupby(cigar_str, lambda c: c.isdigit())
    for _, n in cig_iter:
        yield int("".join(n)), "".join(next(cig_iter)[1])


def read_coverage(cigar_str, read_len):
    """return read coverage as fraction of read length only terminal clips
    and indels are not considered covering

    >>> c="100M"
    >>> read_coverage(c, 100)
    1.0

    >>> c="50S50M"
    >>> read_coverage(c, 100)
    0.5

    >>> c="50M50S"
    >>> read_coverage(c, 100)
    0.5

    >>> c="10S80M10I"
    >>> read_coverage(c, 100)
    0.8

    >>> c="10I80M10S"
    >>> read_coverage(c, 100)
    0.8

    >>> c="5S5I8M80I8M5I5S"
    >>> read_coverage(c, 100)
    0.8

    >>> read_coverage("*", 100)
    0.0
    """

    cigar = list(cigar_items(cigar_str))
    if len(cigar) == 1 and cigar[0][1] is None:
        return 0.0

    covered = read_len
    for (clen, cop) in cigar:
        if cop in "SID":
            covered -= clen
        else:
            break
    for (clen, cop) in reversed(cigar):
        if cop in "SID":
            covered -= clen
        else:
            break
    assert covered >= 0
    return covered/float(read_len)


def min_cov_met(cstr, rlen, mincov):
    """Check whether alignment meets minimum coverage

    >>> min_cov_met("10M", 10, 1.0)
    True
    >>> min_cov_met("5S5M", 10, 1.0)
    False
    >>> min_cov_met("5S5M", 10, 0.5)
    True
    """

    if mincov:
        assert mincov >= 0.0 and mincov <= 1.0
    if not mincov or mincov <= 0.0:
        return True
    if read_coverage(cstr, rlen) < mincov:
        return False
    else:
        return True


def main(fastq_in, ref, fastq_fh, bam_fh, num_threads=2, bwa='bwa', mincov=0.0):
    """main function

    fastq_in and fastq_fh should be lists (single or paired-end)
    """
    # FIXME bufsize needs testing & optimization
    bufsize = 2**16

    assert len(fastq_in) == len(fastq_fh)
    for f in fastq_in + [ref]:
        assert os.path.exists(f)

    bwa_log = tempfile.NamedTemporaryFile(
        mode='w', prefix=os.path.basename(sys.argv[0]) + ".bwa.", delete=False)
    samtools_log = tempfile.NamedTemporaryFile(
        mode='w', prefix=os.path.basename(sys.argv[0]) + ".samtools.", delete=False)

    LOG.info("Using %s as log file for samtools", bwa_log.name)
    LOG.info("Using %s as log file for bwa", samtools_log.name)

    # open BAM as pipe
    # can't use pysam with subprocess
    # p = subprocess.Popen(bwamem_cmd, stdout=subprocess.PIPE)
    # s = pysam.Samfile(p.stdout, "rb")
    # throws 'TypeError: Argument must be string or unicode'
    # see also https://www.biostars.org/p/15298/#105456
    # could use stringio?
    cmd = ["samtools", "view", "-bS", "-"]
    samtools_p = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                  stdout=bam_fh, stderr=samtools_log,
                                  bufsize=bufsize)

    # open BWA as pipe
    #
    # could redirect stderr to file and cat on problem but leave it to
    # user.
    #
    # FIXME add '-M' Mark shorter split hits as secondary. Not for
    # Picard but to make downstream filtering easier (just 0x200
    # instead of 0x200 and 0x800)
    #
    cmd = [bwa, 'mem', '-t', "%d" % num_threads, ref]
    cmd.extend(fastq_in)
    bwa_p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=bwa_log, bufsize=bufsize)

    prev_read = None
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

        counts['Total received from BWA'] += 1
        cur_read = sam_to_read(line)

        # from the SAM spec http://samtools.github.io/hts-specs/SAMv1.pdf:
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

        if cur_read.flag & 0x900:
            counts['Skipped secondary alignments'] += 1
            continue

        is_paired = cur_read.flag & 0x1
        if len(args.fastq_in) == 2:
            assert is_paired, ("Got two fastq files (paired end) but bwa-mem"
                               " output reported single end mapping")
        else:
            assert not is_paired, ("Got one fastq file (single end) but bwa-mem"
                                   " output reported paired end mapping")

        if is_paired and not prev_read:
            # just store for now, no writing
            prev_read = cur_read
            continue

        if is_paired:
            # pairs required to be received sequentially
            assert read_base_name(cur_read.qname) == \
                read_base_name(prev_read.qname), (
                    "Current read name is %s which doesn't match %s" % (
                        cur_read.qname, prev_read.qname))
        else:
            assert not prev_read

        cur_read_mapped = not (cur_read.flag & 0x4)
        if cur_read_mapped and not min_cov_met(
                cur_read.cigar, len(cur_read.seq), mincov):
            cur_read_mapped = False

        reads_to_print = [cur_read]
        if cur_read_mapped:
            to_bam = True
        else:
            to_bam = False


        if is_paired:
            reads_to_print.insert(0, prev_read)

            prev_read_mapped = not (prev_read.flag & 0x4)
            if prev_read_mapped and not min_cov_met(
                    prev_read.cigar, len(prev_read.seq), mincov):
                prev_read_mapped = False

            # one in pair mapped? write both to BAM. oterhwise to fastq
            if cur_read_mapped or prev_read_mapped:
                to_bam = True
            else:
                to_bam = False

            # reset
            prev_read = None


        for r in reads_to_print:
            if to_bam:
                counts['written to %s' % bam_fh.name] += 1
                samtools_p.stdin.write(read_to_sam(r))
            else:
                if not is_paired or r.flag & 0x40:# 1st in pair
                    counts['written to %s' % fastq_fh[0].name] += 1
                    read_to_fastq(r, fastq_fh[0])
                elif r.flag & 0x80:# 2nd in pair
                    counts['written to %s' % fastq_fh[1].name] += 1
                    read_to_fastq(r, fastq_fh[1])
                else:
                    raise ValueError(), ("Read with flag %d neither first nor last in pair" % r.flag)

    for (k, v) in counts.items():
        LOG.info("Reads %s: %d", k, v)

    bwa_p.stdout.close()
    if bwa_p.wait() != 0 or bwa_p.returncode != 0:
        LOG.critical("Unhandled BWA error while processing %s. Check %s",
                     ' and '.join(fastq_in), bwa_log.name)
        return False

    samtools_p.stdin.close()
    if samtools_p.wait() != 0 or samtools_p.returncode != 0:
        LOG.critical("Unhandled samtools error while processing %s. Check %s",
                     ' and '.join(fastq_in), samtools_log.name)
        LOG.critical("Note, this can happen if no reads whatsoever were contaminated.")
        #LOG.critical("samtools might then produce the following complaint:"
        #             " \"reference 'XYZ' is recognized as '*'\""
        #             " followed by \"truncated file\"")
        return False

    return True


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
    mandatory.add_argument("-c", "--mincov",
                           type=float,
                           default=0.0,
                           dest='mincov',
                           help="Minimum alignment coverage to consider a"
                           " read aligned (default: 0.0, i.e. off)")
    default = 8
    parser.add_argument("-t", "--threads",
                        dest='num_threads',
                        default=default,
                        type=int,
                        help="Number of threads to use for mapping"
                        " (default = %d)" % default)
    parser.add_argument("-b", "--bwa",
                        dest='bwa',
                        default="bwa",
                        help="Path to BWA binary (will use BWA found in"
                        " PATH if not set)")
    args = parser.parse_args()

    fastq_out = ["%s_%d.fastq.gz" % (args.outpref, i+1)
                 for i in range(len(args.fastq_in))]
    bam_out = "%s.bam" % (args.outpref)

    infiles = args.fastq_in + [args.ref]
    for f in infiles:
        if not os.path.exists(f):
            LOG.fatal("Input file %s does not exist", f)
            sys.exit(1)

    outfiles = fastq_out + [bam_out]
    for f in outfiles:
        if os.path.exists(f):
            LOG.fatal("Cowardly refusing to overwrite"
                      " already existing file %s", f)
            sys.exit(1)

    fastq_fh = [gzip.open(f, 'w') for f in fastq_out]

    bam_fh = open(bam_out, 'wb')

    if not os.path.exists(args.ref + ".bwt"):
        LOG.warn("Doesn't look like reference was indexed."
                 " Did you forget to run bwa index %s ?", args.ref)

    if not bwa_mem_support(args.bwa):
        LOG.fatal("%s doesn't seem to support mem command.")
        sys.exit(1)

    if args.mincov < 0.0 or args.mincov > 1.0:
        LOG.fatal("minimum coverage arg must be between"
                  " 0 and 1 (but is %f)", args.mincov)
        sys.exit(1)

    rc = main(args.fastq_in, args.ref, fastq_fh, bam_fh,
              num_threads=args.num_threads, bwa=args.bwa,
              mincov=args.mincov)

    for f in fastq_fh + [bam_fh]:
        f.close()

    if not rc:
        sys.exit(1)

    LOG.info("Successful program exit")
