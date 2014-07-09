#!/bin/bash
REF=denv1_refseq.fa
DECONT=../decont.py
OUTPREF=_test
LOG=${OUTPREF}.log
# md5sum is md5 on mac
md5=$(which md5sum 2>/dev/null || which md5)
pushd $(dirname $0)>/dev/null


BWA=bwa
if [ ! -z "$1" ]; then
    BWA=$1
fi

echo "INFO: testing $DECONT with $BWA"
echo "INFO: log file is $LOG and output prefix is always $OUTPREF. If things go wrong, check those files"
echo "INFO: last line should be 'All tests passed successfully'"

# index if needed
test -e ${REF}.bwt || bwa index $REF

# clean up previous unsuccessful run
find . -name ${OUTPREF}\* -exec rm {} \;


test="All contaminated (SR)"
in="SRR1056478_denv1_match_1.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test" 1>&2;
    exit 1;
fi
if [ $n_bam -ne $n_in ] || [ $n_fq -ne 0 ]; then
    echo "ERROR: expected all reads in BAM; none in fastq" 1>&2;
    exit 1;
fi
echo "$test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="All contaminated (PE)"
in="SRR1056478_denv1_match_1.fastq.gz SRR1056478_denv1_match_2.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test" 1>&2;
    exit 1;
fi
if [ $n_bam -ne $n_in ] || [ $n_fq -ne 0 ]; then
    echo "ERROR: expected all reads in BAM; none in fastq" 1>&2;
    exit 1;
fi
echo "$test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="All clean (SR)"
in="SRR1056478_denv1_nomatch_1.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test"; 1>&2
    exit 1;
fi
if [ $n_fq -ne $n_in ] || [ $n_bam -ne 0 ]; then
    echo "ERROR: expected all reads in FastQ; none in BAM" 1>&2;
    exit 1;
fi
md5_in=$(zcat $in | $md5 | cut -f1 -d ' ' )
md5_out=$(zcat ${OUTPREF}_*.fastq.gz | $md5 | cut -f1 -d ' ')
if [ $md5_in != $md5_out ]; then
    echo "ERROR: reads in input and output FastQ differ" 1>&2;
    exit 1;
fi
echo "$test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="All clean (PE)"
in="SRR1056478_denv1_nomatch_1.fastq.gz SRR1056478_denv1_nomatch_2.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test"; 1>&2
    exit 1;
fi
if [ $n_fq -ne $n_in ] || [ $n_bam -ne 0 ]; then
    echo "ERROR: expected all reads in FastQ; none in BAM" 1>&2;
    exit 1;
fi
md5_in=$(zcat $in | $md5 | cut -f1 -d ' ')
md5_out=$(zcat ${OUTPREF}_*.fastq.gz | $md5 | cut -f1 -d ' ')
if [ $md5_in != $md5_out ]; then
    echo "ERROR: reads in input and output FastQ differ" 1>&2;
    exit 1;
fi
echo "$test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


echo "NOTE: missing test for pairs where one maps the other one doesn't" 1>&2

echo "All tests passed successfully"
popd >/dev/null
