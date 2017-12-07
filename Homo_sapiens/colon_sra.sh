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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/colon ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/colon ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030892
mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892

echo "* Start srr_dump [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030892.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` colon ERR030892 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030884
mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884

echo "* Start srr_dump [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030884.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` colon ERR030884 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030892

echo "* Start fastqc [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/ERR030892.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030892 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030884

echo "* Start fastqc [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/ERR030884_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/ERR030884_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030884 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030892

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed

echo "* Start scythe [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/ERR030892.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ERR030892.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ERR030892.scythe.fq.gz


[ $? -ne 0 ] && echo `date` colon ERR030892 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030884

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed

echo "* Start scythe [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/ERR030884_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/ERR030884_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` colon ERR030884 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030892
cd /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed

echo "* Start sickle [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ERR030892.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ERR030892.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030892 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030884
cd /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed

echo "* Start sickle [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030884 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030892

echo "* Start fastqc [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030892/trimmed/ERR030892.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030892 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [colon] [ERR030892] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030884

echo "* Start fastqc [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/colon/ERR030884/trimmed/ERR030884_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` colon ERR030884 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [colon] [ERR030884] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


