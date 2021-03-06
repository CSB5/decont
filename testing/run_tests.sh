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

#  -d no-member for namedtuple FP errors (order matters!)
if ! pylint  -E $DECONT -d no-member; then
        echo "FATAL: pylint failed" 1>&2
        exit 1
fi

if ! python -m doctest $DECONT; then
	echo "FATAL: python doctest failed" 1>&2
	exit 1
fi


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
else
    echo "OK: number of input ($n_in) and output ($n_out) reads match"
fi
if [ $n_bam -ne $n_in ] || [ $n_fq -ne 0 ]; then
    echo "ERROR: expected all reads in BAM; none in fastq" 1>&2;
    exit 1;
else
    echo "OK: all reads in BAM; none in fastq"
fi    
echo "TEST $test: OK"
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
else
    echo "OK: number of input and output reads match"
fi
if [ $n_bam -ne $n_in ] || [ $n_fq -ne 0 ]; then
    echo "ERROR: expected all reads in BAM; none in fastq" 1>&2;
    exit 1;
else
    echo "OK: all reads in BAM; none in fastq"
fi
echo "TEST $test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="All clean (SR)"
in="SRR1056478_denv1_nomatch_1.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
#if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
# "fail" ok, because non contaminated
echo $cmd >$LOG
eval $cmd 2>>$LOG
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test" 1>&2
    exit 1;
else
    echo "OK: number of input and output reads identical"
fi
if [ $n_fq -ne $n_in ] || [ $n_bam -ne 0 ]; then
    echo "ERROR: expected all reads in FastQ; none in BAM" 1>&2;
    exit 1;
else
    echo "OK: all reads in FastQ; none in BAM"
fi
md5_in=$(gzip -dc $in | $md5 | cut -f1 -d ' ' )
md5_out=$(gzip -dc ${OUTPREF}_*.fastq.gz | $md5 | cut -f1 -d ' ')
if [ $md5_in != $md5_out ]; then
    echo "ERROR: reads in input ($in) and output FastQ $(ls ${OUTPREF}_*.fastq.gz) differ" 1>&2;
    exit 1;
else
    echo "OK: reads in input ($in) and output FastQ identical"
fi
    
echo "TEST $test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;




test="All clean (PE)"
in="SRR1056478_denv1_nomatch_1.fastq.gz SRR1056478_denv1_nomatch_2.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
#if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
# "fail" ok, because non contaminated
eval $cmd 2>$LOG
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
n_out=$(echo $n_bam+$n_fq | bc)
#echo "DEBUG: $test: n_in=$n_in n_bam=$n_bam n_fq=$n_fq n_out=$n_out" 1>&2
if [ $n_in -ne $n_out ]; then
    echo "ERROR: mismatch between number of input ($n_in) and output ($n_out) reads in test $test" 1>&2
    exit 1;
else
    echo "OK: number of input and output reads match"
fi
if [ $n_fq -ne $n_in ] || [ $n_bam -ne 0 ]; then
    echo "ERROR: expected all reads in FastQ; none in BAM" 1>&2;
    exit 1;
else
    echo "OK: all reads in FastQ; none in BAM"
fi 
md5_in=$(gzip -dc $in | $md5 | cut -f1 -d ' ')
md5_out=$(gzip -dc ${OUTPREF}_*.fastq.gz | $md5 | cut -f1 -d ' ')
if [ $md5_in != $md5_out ]; then
    echo "ERROR: reads in input ($in) and output $(ls ${OUTPREF}_*.fastq.gz) FastQ differ" 1>&2;
    exit 1;
else
    echo "OK: reads in input and output FastQ identical"
fi  
echo "TEST $test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="Match/Nonmatch pair"
in="pe_match.fastq.gz pe_nomatch.fastq.gz"
cmd="$DECONT -b $BWA -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
if [ $n_in -ne $n_bam ]; then
    echo "ERROR: expected match/non-match pair to end up in BAM but didn't in $test" 1>&2
    exit 1;
else
    echo "OK: match/non-match pair ended up in BAM"
fi
echo "TEST $test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="Minimum coverage"
in="pe_match_50cov.fastq.gz"
cmd="$DECONT -b $BWA -c 0.51 -i $in -o ${OUTPREF} -r $REF"
#if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
# "fail" ok, because non contaminated
eval $cmd 2>$LOG
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_fq=$(fastq_num_reads.sh ${OUTPREF}_*.fastq.gz | awk '{s+=$NF} END {print s}')
if [ $n_in -ne $n_fq ]; then
    echo "ERROR: expected matches below coverage limit to count as unaligned" 1>&2
    exit 1;
else
    echo "OK matches below coverage limit counted as unaligned"
fi
find . -name ${OUTPREF}\* -exec rm {} \;
cmd="$DECONT -b $BWA -c 0.49 -i $in -o ${OUTPREF} -r $REF"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
n_in=$(fastq_num_reads.sh $in | awk '{s+=$NF} END {print s}')
n_bam=$(samtools view -c ${OUTPREF}.bam)
if [ $n_in -ne $n_bam ]; then
    echo "ERROR: expected matches above coverage limit to count as aligned" 1>&2
    exit 1;
else
    echo "OK: matches above coverage limit counted as aligned"
fi
echo "TEST $test: OK"
find . -name ${OUTPREF}\* -exec rm {} \;


test="Catch errors"
# if bwa eats it, we will
#in="invalid.fastq"
#cmd="$DECONT -i $in -r $REF -o ${OUTPREF}"
#if eval $cmd 2>$LOG; then echo "ERROR: the following command should have failed, but didn't: $cmd" 1>&2; exit 1; fi
#echo "$test: OK"
#find . -name ${OUTPREF}\* -exec rm {} \;


test="Pair-suffix present in FastQ, but not in BAM"
in="SRR1056478_denv1_match_1.fastq.gz SRR1056478_denv1_match_2.fastq.gz"
cmd="$DECONT -i $in -r $REF -o ${OUTPREF}"
if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
name=$(samtools view ${OUTPREF}.bam | head -n 1 | cut -f 1) || exit 1
if echo $name | grep -q '[#/][12]$'; then
    echo "ERROR: pair-suffices should be missing in BAM but aren't" 1>&2
else
    echo "OK: pair-suffices are missing in BAM"
fi
find . -name ${OUTPREF}\* -exec rm {} \;
in="SRR1056478_denv1_nomatch_1.fastq.gz SRR1056478_denv1_nomatch_2.fastq.gz"
cmd="$DECONT -i $in -r $REF -o ${OUTPREF}"
#if ! eval $cmd 2>$LOG; then echo "ERROR: the following command failed: $cmd" 1>&2; exit 1; fi
# "fail" ok, because non contaminated
eval $cmd 2>$LOG
name=$(gzip -dc $in | head -n 1)
if ! echo $name | grep -q '[#/][12]$'; then
    echo "ERROR: pair-suffices should be present in FastQ but aren't" 1>&2
else
    echo "OK: pair-suffices are present in FastQ"
fi
find . -name ${OUTPREF}\* -exec rm {} \;


echo "All tests passed successfully"
find . -name ${OUTPREF}\* -exec rm {} \;
popd >/dev/null
