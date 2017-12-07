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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/testes ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/testes ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030902
mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902

echo "* Start srr_dump [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030902.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` testes ERR030902 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030873
mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873

echo "* Start srr_dump [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030873.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` testes ERR030873 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030902

echo "* Start fastqc [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/ERR030902.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030902 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030873

echo "* Start fastqc [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/ERR030873_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/ERR030873_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030873 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030902

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed

echo "* Start scythe [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/ERR030902.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ERR030902.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ERR030902.scythe.fq.gz


[ $? -ne 0 ] && echo `date` testes ERR030902 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030873

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed

echo "* Start scythe [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/ERR030873_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/ERR030873_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` testes ERR030873 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030902
cd /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed

echo "* Start sickle [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ERR030902.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ERR030902.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030902 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030873
cd /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed

echo "* Start sickle [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030873 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030902

echo "* Start fastqc [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030902/trimmed/ERR030902.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030902 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [testes] [ERR030902] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030873

echo "* Start fastqc [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/testes/ERR030873/trimmed/ERR030873_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` testes ERR030873 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [testes] [ERR030873] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


