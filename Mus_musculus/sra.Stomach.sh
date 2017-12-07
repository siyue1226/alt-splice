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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Stomach ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453093
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093

echo "* Start srr_dump [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453093.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453093 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453094
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094

echo "* Start srr_dump [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453094.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453094 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453095
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095

echo "* Start srr_dump [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453095.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453095 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453096
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096

echo "* Start srr_dump [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453096.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453096 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453097
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097

echo "* Start srr_dump [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453097.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453097 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453098
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098

echo "* Start srr_dump [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453098.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Stomach SRR453098 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453093

echo "* Start fastqc [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/SRR453093_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/SRR453093_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453093 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453094

echo "* Start fastqc [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/SRR453094_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/SRR453094_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453094 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453095

echo "* Start fastqc [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/SRR453095_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/SRR453095_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453095 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453096

echo "* Start fastqc [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/SRR453096_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/SRR453096_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453096 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453097

echo "* Start fastqc [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/SRR453097_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/SRR453097_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453097 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453098

echo "* Start fastqc [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/SRR453098_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/SRR453098_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453098 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453093

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed

echo "* Start scythe [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/SRR453093_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/SRR453093_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453093 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453094

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed

echo "* Start scythe [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/SRR453094_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/SRR453094_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453094 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453095

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed

echo "* Start scythe [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/SRR453095_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/SRR453095_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453095 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453096

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed

echo "* Start scythe [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/SRR453096_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/SRR453096_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453096 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453097

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed

echo "* Start scythe [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/SRR453097_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/SRR453097_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453097 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453098

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed

echo "* Start scythe [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/SRR453098_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/SRR453098_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453098 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453093
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed

echo "* Start sickle [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453093 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453094
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed

echo "* Start sickle [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453094 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453095
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed

echo "* Start sickle [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453095 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453096
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed

echo "* Start sickle [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453096 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453097
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed

echo "* Start sickle [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453097 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453098
cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed

echo "* Start sickle [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453098 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453093

echo "* Start fastqc [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453093 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453094

echo "* Start fastqc [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453094 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453095

echo "* Start fastqc [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453095 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453096

echo "* Start fastqc [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453096 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453097

echo "* Start fastqc [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453097 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453098

echo "* Start fastqc [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Stomach SRR453098 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


