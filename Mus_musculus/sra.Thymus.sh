#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/mouse_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa
# samtools faidx /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Thymus ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453134
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134

echo "* Start srr_dump [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453134.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453134 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453135
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135

echo "* Start srr_dump [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453135.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453135 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453136
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136

echo "* Start srr_dump [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453136.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453136 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453137
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137

echo "* Start srr_dump [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453137.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453137 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453138
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138

echo "* Start srr_dump [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453138.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453138 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453139
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139

echo "* Start srr_dump [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453139.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Thymus SRR453139 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453134

echo "* Start fastqc [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/SRR453134_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/SRR453134_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453134 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453135

echo "* Start fastqc [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/SRR453135_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/SRR453135_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453135 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453136

echo "* Start fastqc [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/SRR453136_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/SRR453136_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453136 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453137

echo "* Start fastqc [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/SRR453137_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/SRR453137_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453137 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453138

echo "* Start fastqc [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/SRR453138_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/SRR453138_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453138 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453139

echo "* Start fastqc [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/SRR453139_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/SRR453139_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453139 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453134

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed

echo "* Start scythe [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/SRR453134_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/SRR453134_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453134 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453135

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed

echo "* Start scythe [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/SRR453135_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/SRR453135_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453135 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453136

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed

echo "* Start scythe [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/SRR453136_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/SRR453136_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453136 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453137

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed

echo "* Start scythe [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/SRR453137_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/SRR453137_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453137 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453138

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed

echo "* Start scythe [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/SRR453138_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/SRR453138_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453138 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453139

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed

echo "* Start scythe [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/SRR453139_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/SRR453139_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453139 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453134
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed

echo "* Start sickle [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453134 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453135
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed

echo "* Start sickle [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453135 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453136
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed

echo "* Start sickle [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453136 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453137
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed

echo "* Start sickle [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453137 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453138
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed

echo "* Start sickle [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453138 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453139
cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed

echo "* Start sickle [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453139 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453134

echo "* Start fastqc [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453134 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453135

echo "* Start fastqc [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453135 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453136

echo "* Start fastqc [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453136 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453137

echo "* Start fastqc [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453137 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453138

echo "* Start fastqc [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453138 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453139

echo "* Start fastqc [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Thymus SRR453139 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


