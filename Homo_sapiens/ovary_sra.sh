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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/ovary ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/ovary ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030874
mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874

echo "* Start srr_dump [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030874.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` ovary ERR030874 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030901
mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901

echo "* Start srr_dump [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030901.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` ovary ERR030901 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030874

echo "* Start fastqc [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/ERR030874_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/ERR030874_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030874 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030901

echo "* Start fastqc [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/ERR030901.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030901 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030874

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed

echo "* Start scythe [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/ERR030874_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/ERR030874_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` ovary ERR030874 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030901

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed

echo "* Start scythe [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/ERR030901.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.scythe.fq.gz


[ $? -ne 0 ] && echo `date` ovary ERR030901 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030874
cd /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed

echo "* Start sickle [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030874 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030901
cd /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed

echo "* Start sickle [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030901 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030874

echo "* Start fastqc [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030874 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [ovary] [ERR030874] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030901

echo "* Start fastqc [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` ovary ERR030901 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [ovary] [ERR030901] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


