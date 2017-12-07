#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq1/data/rna2dna/bodymap2

### index reference genome
# bwa index -a bwtsw /home/wangq1/data/rna2dna/bodymap2/ref/human.37.fa
# samtools faidx /home/wangq1/data/rna2dna/bodymap2/ref/human.37.fa

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/brain ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/brain ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030890
mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890

echo "* Start srr_dump [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030890.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` brain ERR030890 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030882
mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882

echo "* Start srr_dump [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030882.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` brain ERR030882 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030890

echo "* Start fastqc [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/ERR030890.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030890 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030882

echo "* Start fastqc [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/ERR030882_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/ERR030882_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030882 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030890

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed

echo "* Start scythe [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/ERR030890.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ERR030890.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ERR030890.scythe.fq.gz


[ $? -ne 0 ] && echo `date` brain ERR030890 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030882

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed

echo "* Start scythe [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/ERR030882_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/ERR030882_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` brain ERR030882 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030890
cd /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed

echo "* Start sickle [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ERR030890.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ERR030890.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030890 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030882
cd /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed

echo "* Start sickle [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030882 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030890

echo "* Start fastqc [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030890/trimmed/ERR030890.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030890 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [brain] [ERR030890] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030882

echo "* Start fastqc [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/brain/ERR030882/trimmed/ERR030882_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` brain ERR030882 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [brain] [ERR030882] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


