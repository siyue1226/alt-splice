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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/kidney ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/kidney ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030885
mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885

echo "* Start srr_dump [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030885.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` kidney ERR030885 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030893
mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893

echo "* Start srr_dump [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030893.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` kidney ERR030893 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030885

echo "* Start fastqc [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/ERR030885_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/ERR030885_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030885 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030893

echo "* Start fastqc [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/ERR030893.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030893 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030885

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed

echo "* Start scythe [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/ERR030885_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/ERR030885_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` kidney ERR030885 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030893

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed

echo "* Start scythe [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/ERR030893.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.scythe.fq.gz


[ $? -ne 0 ] && echo `date` kidney ERR030893 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030885
cd /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed

echo "* Start sickle [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030885 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030893
cd /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed

echo "* Start sickle [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030893 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030885

echo "* Start fastqc [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030885 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [kidney] [ERR030885] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030893

echo "* Start fastqc [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` kidney ERR030893 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [kidney] [ERR030893] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


