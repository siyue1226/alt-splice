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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/lung ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/lung ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030896
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896

echo "* Start srr_dump [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030896.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` lung ERR030896 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030879
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879

echo "* Start srr_dump [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030879.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` lung ERR030879 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030896

echo "* Start fastqc [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/ERR030896.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030896 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030879

echo "* Start fastqc [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/ERR030879_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/ERR030879_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030879 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030896

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed

echo "* Start scythe [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/ERR030896.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ERR030896.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ERR030896.scythe.fq.gz


[ $? -ne 0 ] && echo `date` lung ERR030896 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030879

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed

echo "* Start scythe [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/ERR030879_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/ERR030879_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` lung ERR030879 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030896
cd /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed

echo "* Start sickle [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ERR030896.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ERR030896.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030896 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030879
cd /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed

echo "* Start sickle [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030879 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030896

echo "* Start fastqc [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030896/trimmed/ERR030896.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030896 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lung] [ERR030896] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030879

echo "* Start fastqc [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lung/ERR030879/trimmed/ERR030879_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lung ERR030879 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lung] [ERR030879] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


