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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/liver ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/liver ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030895
mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895

echo "* Start srr_dump [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030895.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` liver ERR030895 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030887
mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887

echo "* Start srr_dump [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030887.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` liver ERR030887 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030895

echo "* Start fastqc [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/ERR030895.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030895 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030887

echo "* Start fastqc [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/ERR030887_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/ERR030887_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030887 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030895

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed

echo "* Start scythe [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/ERR030895.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ERR030895.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ERR030895.scythe.fq.gz


[ $? -ne 0 ] && echo `date` liver ERR030895 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030887

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed

echo "* Start scythe [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/ERR030887_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/ERR030887_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` liver ERR030887 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030895
cd /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed

echo "* Start sickle [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ERR030895.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ERR030895.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030895 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030887
cd /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed

echo "* Start sickle [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030887 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030895

echo "* Start fastqc [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030895/trimmed/ERR030895.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030895 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [liver] [ERR030895] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030887

echo "* Start fastqc [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/liver/ERR030887/trimmed/ERR030887_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` liver ERR030887 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [liver] [ERR030887] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


