#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/rice_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa
# samtools faidx /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa

if [ -d /home/wangq/data/rna-seq/rice_trans/process/ ];
then
    rm -fr /home/wangq/data/rna-seq/rice_trans/process/ ;
fi;
mkdir /home/wangq/data/rna-seq/rice_trans/process/

#----------------------------#
# srr dump
#----------------------------#
# lane SRR352207
mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352207

echo "* Start srr_dump [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/rice_trans/SRR352207.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/rice_trans/process/SRR352207 \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` OSN_AK SRR352207 [fastq dump] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End srr_dump [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR352207

echo "* Start fastqc [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/SRR352207_1.fastq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/SRR352207_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AK SRR352207 [fastqc] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End fastqc [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log



#----------------------------#
# scythe
#----------------------------#
# lane SRR352207

if [ ! -d /home/wangq/data/rna-seq/rice_trans/process/SRR352207 ];
then
    mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352207;
fi;

if [ ! -d /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed ;
fi;

cd /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed

echo "* Start scythe [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/SRR352207_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/rice_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/SRR352207_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/rice_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` OSN_AK SRR352207 [scythe] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End scythe [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/scythe.log



#----------------------------#
# sickle
#----------------------------#
# lane SRR352207
cd /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed

echo "* Start sickle [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_1.sickle.fq \
    -p /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_2.sickle.fq \
    -s /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AK SRR352207 [sickle] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End sickle [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log

find /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR352207

echo "* Start fastqc [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/trimmed/SRR352207_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AK SRR352207 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End fastqc [OSN_AK] [SRR352207] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log



