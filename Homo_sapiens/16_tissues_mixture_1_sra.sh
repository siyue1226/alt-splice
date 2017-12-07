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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1 ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1 ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030868
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030868.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030869
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030869.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030871
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030871.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030863
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030863.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030862
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030862.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030870
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870

echo "* Start srr_dump [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030870.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030868

echo "* Start fastqc [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/ERR030868.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030869

echo "* Start fastqc [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/ERR030869.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030871

echo "* Start fastqc [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/ERR030871.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030863

echo "* Start fastqc [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/ERR030863.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030862

echo "* Start fastqc [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/ERR030862.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030870

echo "* Start fastqc [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/ERR030870.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030868

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/ERR030868.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030869

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/ERR030869.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030871

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/ERR030871.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030863

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/ERR030863.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030862

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/ERR030862.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030870

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed

echo "* Start scythe [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/ERR030870.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030868
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030869
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030871
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030863
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030862
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030870
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed

echo "* Start sickle [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030868

echo "* Start fastqc [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030869

echo "* Start fastqc [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030871

echo "* Start fastqc [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030863

echo "* Start fastqc [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030862

echo "* Start fastqc [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030870

echo "* Start fastqc [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


