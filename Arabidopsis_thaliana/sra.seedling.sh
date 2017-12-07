#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/Ath_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa
# samtools faidx /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa

if [ -d /home/wangq/data/rna-seq/Ath_trans/process/ ];
then
    rm -fr /home/wangq/data/rna-seq/Ath_trans/process/ ;
fi;
mkdir /home/wangq/data/rna-seq/Ath_trans/process/

#----------------------------#
# srr dump
#----------------------------#
# lane SRR6179908
mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908

echo "* Start srr_dump [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/Ath_trans/SRR6179908.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908 \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` seedling SRR6179908 [fastq dump] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End srr_dump [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR6179908

echo "* Start fastqc [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/SRR6179908_1.fastq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/SRR6179908_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` seedling SRR6179908 [fastqc] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End fastqc [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log



#----------------------------#
# scythe
#----------------------------#
# lane SRR6179908

if [ ! -d /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908 ];
then
    mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908;
fi;

if [ ! -d /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed ;
fi;

cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed

echo "* Start scythe [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/SRR6179908_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/Ath_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/SRR6179908_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/Ath_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` seedling SRR6179908 [scythe] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End scythe [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/scythe.log



#----------------------------#
# sickle
#----------------------------#
# lane SRR6179908
cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed

echo "* Start sickle [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_1.sickle.fq \
    -p /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_2.sickle.fq \
    -s /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` seedling SRR6179908 [sickle] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End sickle [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log

find /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR6179908

echo "* Start fastqc [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/trimmed/SRR6179908_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` seedling SRR6179908 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End fastqc [seedling] [SRR6179908] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log



