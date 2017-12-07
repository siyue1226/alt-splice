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
# lane SRR6179905
mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905

echo "* Start srr_dump [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/Ath_trans/SRR6179905.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905 \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` steam SRR6179905 [fastq dump] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End srr_dump [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/srr_dump.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR6179905

echo "* Start fastqc [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/SRR6179905_1.fastq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/SRR6179905_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` steam SRR6179905 [fastqc] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End fastqc [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log



#----------------------------#
# scythe
#----------------------------#
# lane SRR6179905

if [ ! -d /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905 ];
then
    mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905;
fi;

if [ ! -d /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed ;
fi;

cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed

echo "* Start scythe [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/SRR6179905_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/Ath_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/SRR6179905_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/Ath_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` steam SRR6179905 [scythe] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End scythe [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/scythe.log



#----------------------------#
# sickle
#----------------------------#
# lane SRR6179905
cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed

echo "* Start sickle [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_1.sickle.fq \
    -p /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_2.sickle.fq \
    -s /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` steam SRR6179905 [sickle] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End sickle [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log

find /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/sickle.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR6179905

echo "* Start fastqc [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/trimmed/SRR6179905_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` steam SRR6179905 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End fastqc [steam] [SRR6179905] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/fastqc.log



