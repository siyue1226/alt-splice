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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030875
mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875

echo "* Start srr_dump [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030875.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030900
mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900

echo "* Start srr_dump [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030900.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030875

echo "* Start fastqc [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/ERR030875_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/ERR030875_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030900

echo "* Start fastqc [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/ERR030900.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030875

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed

echo "* Start scythe [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/ERR030875_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/ERR030875_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030900

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed

echo "* Start scythe [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/ERR030900.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.scythe.fq.gz


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030875
cd /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed

echo "* Start sickle [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030900
cd /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed

echo "* Start sickle [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030875

echo "* Start fastqc [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030900

echo "* Start fastqc [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


