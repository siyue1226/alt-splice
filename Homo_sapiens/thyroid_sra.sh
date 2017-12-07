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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/thyroid ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/thyroid ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030872
mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872

echo "* Start srr_dump [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030872.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` thyroid ERR030872 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030903
mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903

echo "* Start srr_dump [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030903.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` thyroid ERR030903 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030872

echo "* Start fastqc [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/ERR030872_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/ERR030872_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030872 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030903

echo "* Start fastqc [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/ERR030903.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030903 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030872

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed

echo "* Start scythe [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/ERR030872_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/ERR030872_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` thyroid ERR030872 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030903

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed

echo "* Start scythe [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/ERR030903.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.scythe.fq.gz


[ $? -ne 0 ] && echo `date` thyroid ERR030903 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030872
cd /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed

echo "* Start sickle [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030872 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030903
cd /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed

echo "* Start sickle [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030903 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030872

echo "* Start fastqc [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030872 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [thyroid] [ERR030872] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030903

echo "* Start fastqc [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` thyroid ERR030903 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [thyroid] [ERR030903] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


