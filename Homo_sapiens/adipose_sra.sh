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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/adipose ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/adipose ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030888
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888

echo "* Start srr_dump [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030888.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` adipose ERR030888 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030880
mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880

echo "* Start srr_dump [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030880.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` adipose ERR030880 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030888

echo "* Start fastqc [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/ERR030888.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030888 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030880

echo "* Start fastqc [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/ERR030880_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/ERR030880_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030880 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030888

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed

echo "* Start scythe [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/ERR030888.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.scythe.fq.gz


[ $? -ne 0 ] && echo `date` adipose ERR030888 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030880

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed

echo "* Start scythe [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/ERR030880_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/ERR030880_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` adipose ERR030880 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030888
cd /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed

echo "* Start sickle [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030888 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030880
cd /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed

echo "* Start sickle [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030880 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030888

echo "* Start fastqc [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030888 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adipose] [ERR030888] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030880

echo "* Start fastqc [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` adipose ERR030880 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [adipose] [ERR030880] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


