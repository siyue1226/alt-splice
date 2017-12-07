#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/dmel_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa
# samtools faidx /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070422
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422

echo "* Start srr_dump [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070422.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070423
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423

echo "* Start srr_dump [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070423.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100276
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276

echo "* Start srr_dump [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100276.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR350960
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960

echo "* Start srr_dump [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR350960.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR350961
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961

echo "* Start srr_dump [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR350961.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070422

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/SRR070422_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/SRR070422_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070423

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/SRR070423_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/SRR070423_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100276

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/SRR100276_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/SRR100276_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350960

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/SRR350960_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/SRR350960_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350961

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/SRR350961_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/SRR350961_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070422

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed

echo "* Start scythe [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/SRR070422_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/SRR070422_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070423

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed

echo "* Start scythe [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/SRR070423_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/SRR070423_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100276

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed

echo "* Start scythe [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/SRR100276_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/SRR100276_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR350960

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed

echo "* Start scythe [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/SRR350960_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/SRR350960_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR350961

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed

echo "* Start scythe [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/SRR350961_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/SRR350961_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070422
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed

echo "* Start sickle [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070423
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed

echo "* Start sickle [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100276
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed

echo "* Start sickle [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR350960
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed

echo "* Start sickle [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR350961
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed

echo "* Start sickle [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070422

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070423

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100276

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350960

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350961

echo "* Start fastqc [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


