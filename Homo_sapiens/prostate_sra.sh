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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/prostate ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/prostate ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030877
mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877

echo "* Start srr_dump [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030877.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` prostate ERR030877 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030898
mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898

echo "* Start srr_dump [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030898.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` prostate ERR030898 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030877

echo "* Start fastqc [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/ERR030877_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/ERR030877_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030877 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030898

echo "* Start fastqc [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/ERR030898.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030898 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030877

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed

echo "* Start scythe [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/ERR030877_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/ERR030877_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` prostate ERR030877 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030898

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed

echo "* Start scythe [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/ERR030898.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.scythe.fq.gz


[ $? -ne 0 ] && echo `date` prostate ERR030898 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030877
cd /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed

echo "* Start sickle [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030877 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030898
cd /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed

echo "* Start sickle [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030898 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030877

echo "* Start fastqc [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030877 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [prostate] [ERR030877] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030898

echo "* Start fastqc [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` prostate ERR030898 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [prostate] [ERR030898] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


