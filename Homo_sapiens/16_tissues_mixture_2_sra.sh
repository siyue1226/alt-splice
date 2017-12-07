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

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2 ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2 ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030859
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859

echo "* Start srr_dump [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030859.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030867
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867

echo "* Start srr_dump [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030867.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030861
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861

echo "* Start srr_dump [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030861.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030866
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866

echo "* Start srr_dump [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030866.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030860
mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860

echo "* Start srr_dump [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030860.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030859

echo "* Start fastqc [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/ERR030859.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030867

echo "* Start fastqc [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/ERR030867.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030861

echo "* Start fastqc [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/ERR030861.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030866

echo "* Start fastqc [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/ERR030866.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030860

echo "* Start fastqc [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/ERR030860.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030859

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed

echo "* Start scythe [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/ERR030859.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030867

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed

echo "* Start scythe [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/ERR030867.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030861

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed

echo "* Start scythe [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/ERR030861.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030866

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed

echo "* Start scythe [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/ERR030866.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030860

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed

echo "* Start scythe [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/ERR030860.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.scythe.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030859
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed

echo "* Start sickle [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030867
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed

echo "* Start sickle [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030861
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed

echo "* Start sickle [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030866
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed

echo "* Start sickle [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030860
cd /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed

echo "* Start sickle [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030859

echo "* Start fastqc [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030867

echo "* Start fastqc [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030861

echo "* Start fastqc [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030866

echo "* Start fastqc [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030860

echo "* Start fastqc [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


