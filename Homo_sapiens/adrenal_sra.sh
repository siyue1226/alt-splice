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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/adrenal ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/adrenal ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030889
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889

echo "* Start srr_dump [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030889.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` adrenal ERR030889 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030881
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881

echo "* Start srr_dump [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030881.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` adrenal ERR030881 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030889

echo "* Start fastqc [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/ERR030889.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030889 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030881

echo "* Start fastqc [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/ERR030881_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/ERR030881_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030881 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030889

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed

echo "* Start scythe [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/ERR030889.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ERR030889.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ERR030889.scythe.fq.gz


[ $? -ne 0 ] && echo `date` adrenal ERR030889 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030881

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed

echo "* Start scythe [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/ERR030881_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/ERR030881_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` adrenal ERR030881 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030889
cd /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed

echo "* Start sickle [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ERR030889.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ERR030889.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030889 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030881
cd /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed

echo "* Start sickle [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030881 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030889

echo "* Start fastqc [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030889/trimmed/ERR030889.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030889 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adrenal] [ERR030889] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030881

echo "* Start fastqc [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adrenal/ERR030881/trimmed/ERR030881_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adrenal ERR030881 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adrenal] [ERR030881] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


