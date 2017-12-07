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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Adrenal ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453116
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116

echo "* Start srr_dump [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453116.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453117
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117

echo "* Start srr_dump [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453117.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453118
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118

echo "* Start srr_dump [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453118.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453119
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119

echo "* Start srr_dump [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453119.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453120
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120

echo "* Start srr_dump [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453120.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453121
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121

echo "* Start srr_dump [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453121.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453116

echo "* Start fastqc [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/SRR453116_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/SRR453116_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453117

echo "* Start fastqc [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/SRR453117_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/SRR453117_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453118

echo "* Start fastqc [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/SRR453118_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/SRR453118_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453119

echo "* Start fastqc [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/SRR453119_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/SRR453119_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453120

echo "* Start fastqc [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/SRR453120_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/SRR453120_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453121

echo "* Start fastqc [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/SRR453121_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/SRR453121_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453116

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed

echo "* Start scythe [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/SRR453116_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/SRR453116_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453117

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed

echo "* Start scythe [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/SRR453117_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/SRR453117_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453118

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed

echo "* Start scythe [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/SRR453118_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/SRR453118_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453119

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed

echo "* Start scythe [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/SRR453119_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/SRR453119_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453120

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed

echo "* Start scythe [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/SRR453120_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/SRR453120_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453121

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed

echo "* Start scythe [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/SRR453121_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/SRR453121_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453116
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed

echo "* Start sickle [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453117
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed

echo "* Start sickle [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453118
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed

echo "* Start sickle [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453119
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed

echo "* Start sickle [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453120
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed

echo "* Start sickle [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453121
cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed

echo "* Start sickle [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453116

echo "* Start fastqc [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453117

echo "* Start fastqc [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453118

echo "* Start fastqc [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453119

echo "* Start fastqc [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453120

echo "* Start fastqc [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453121

echo "* Start fastqc [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


