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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/breast ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/breast ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030891
mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891

echo "* Start srr_dump [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030891.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` breast ERR030891 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030883
mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883

echo "* Start srr_dump [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030883.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` breast ERR030883 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030891

echo "* Start fastqc [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/ERR030891.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030891 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030883

echo "* Start fastqc [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/ERR030883_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/ERR030883_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030883 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030891

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed

echo "* Start scythe [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/ERR030891.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ERR030891.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ERR030891.scythe.fq.gz


[ $? -ne 0 ] && echo `date` breast ERR030891 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030883

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed

echo "* Start scythe [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/ERR030883_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/ERR030883_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` breast ERR030883 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030891
cd /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed

echo "* Start sickle [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ERR030891.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ERR030891.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030891 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030883
cd /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed

echo "* Start sickle [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030883 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030891

echo "* Start fastqc [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030891/trimmed/ERR030891.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030891 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [breast] [ERR030891] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030883

echo "* Start fastqc [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/breast/ERR030883/trimmed/ERR030883_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` breast ERR030883 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [breast] [ERR030883] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


