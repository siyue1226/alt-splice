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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/heart ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/heart ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030894
mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894

echo "* Start srr_dump [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030894.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` heart ERR030894 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030886
mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886

echo "* Start srr_dump [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030886.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` heart ERR030886 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030894

echo "* Start fastqc [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/ERR030894.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030894 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030886

echo "* Start fastqc [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/ERR030886_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/ERR030886_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030886 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030894

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed

echo "* Start scythe [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/ERR030894.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ERR030894.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ERR030894.scythe.fq.gz


[ $? -ne 0 ] && echo `date` heart ERR030894 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030886

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed

echo "* Start scythe [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/ERR030886_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/ERR030886_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` heart ERR030886 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030894
cd /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed

echo "* Start sickle [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ERR030894.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ERR030894.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030894 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030886
cd /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed

echo "* Start sickle [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030886 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030894

echo "* Start fastqc [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030894/trimmed/ERR030894.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030894 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [heart] [ERR030894] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030886

echo "* Start fastqc [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/heart/ERR030886/trimmed/ERR030886_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` heart ERR030886 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [heart] [ERR030886] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


