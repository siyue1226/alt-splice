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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3 ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3 ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030864
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864

echo "* Start srr_dump [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030864.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030865
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865

echo "* Start srr_dump [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030865.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030857
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857

echo "* Start srr_dump [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030857.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030858
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858

echo "* Start srr_dump [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030858.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030856
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856

echo "* Start srr_dump [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030856.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030864

echo "* Start fastqc [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/ERR030864.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030865

echo "* Start fastqc [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/ERR030865.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030857

echo "* Start fastqc [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/ERR030857.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030858

echo "* Start fastqc [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/ERR030858.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030856

echo "* Start fastqc [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/ERR030856.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030864

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed

echo "* Start scythe [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/ERR030864.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030865

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed

echo "* Start scythe [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/ERR030865.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030857

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed

echo "* Start scythe [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/ERR030857.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030858

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed

echo "* Start scythe [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/ERR030858.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030856

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed

echo "* Start scythe [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/ERR030856.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030864
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed

echo "* Start sickle [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030865
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed

echo "* Start sickle [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030857
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed

echo "* Start sickle [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030858
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed

echo "* Start sickle [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030856
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed

echo "* Start sickle [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030864

echo "* Start fastqc [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030865

echo "* Start fastqc [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030857

echo "* Start fastqc [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030858

echo "* Start fastqc [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030856

echo "* Start fastqc [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


