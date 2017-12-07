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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030876
mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876

echo "* Start srr_dump [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030876.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030899
mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899

echo "* Start srr_dump [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030899.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030876

echo "* Start fastqc [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/ERR030876_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/ERR030876_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030899

echo "* Start fastqc [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/ERR030899.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030876

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed

echo "* Start scythe [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/ERR030876_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/ERR030876_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030899

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed

echo "* Start scythe [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/ERR030899.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.scythe.fq.gz


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030876
cd /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed

echo "* Start sickle [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030899
cd /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed

echo "* Start sickle [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030876

echo "* Start fastqc [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030899

echo "* Start fastqc [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


