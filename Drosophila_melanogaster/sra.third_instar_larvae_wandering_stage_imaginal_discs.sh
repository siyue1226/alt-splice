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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070392
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070392.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070393
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070393.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111884
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111884.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111885
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111885.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR350962
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR350962.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR350963
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963

echo "* Start srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR350963.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070392

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/SRR070392_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/SRR070392_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070393

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/SRR070393_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/SRR070393_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111884

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/SRR111884_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/SRR111884_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111885

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/SRR111885_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/SRR111885_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350962

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/SRR350962_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/SRR350962_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350963

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/SRR350963_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/SRR350963_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070392

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/SRR070392_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/SRR070392_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070393

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/SRR070393_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/SRR070393_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111884

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/SRR111884_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/SRR111884_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111885

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/SRR111885_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/SRR111885_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR350962

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/SRR350962_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/SRR350962_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR350963

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/SRR350963_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/SRR350963_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070392
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070393
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111884
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111885
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR350962
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR350963
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070392

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070393

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111884

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111885

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350962

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350963

echo "* Start fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


