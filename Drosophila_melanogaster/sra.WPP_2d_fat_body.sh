#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/dmel_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa
# samtools faidx /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070429
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429

echo "* Start srr_dump [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070429.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` WPP_2d_fat_body SRR070429 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070429

echo "* Start fastqc [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/SRR070429_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/SRR070429_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_fat_body SRR070429 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070429

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed

echo "* Start scythe [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/SRR070429_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/SRR070429_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` WPP_2d_fat_body SRR070429 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070429
cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed

echo "* Start sickle [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_fat_body SRR070429 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070429

echo "* Start fastqc [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/trimmed/SRR070429_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_fat_body SRR070429 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_fat_body] [SRR070429] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


