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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Duodenum ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453109
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109

echo "* Start srr_dump [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453109.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453110
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110

echo "* Start srr_dump [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453110.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453111
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111

echo "* Start srr_dump [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453111.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453112
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112

echo "* Start srr_dump [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453112.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453113
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113

echo "* Start srr_dump [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453113.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453114
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114

echo "* Start srr_dump [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453114.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453115
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115

echo "* Start srr_dump [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453115.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453109

echo "* Start fastqc [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/SRR453109_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/SRR453109_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453110

echo "* Start fastqc [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/SRR453110_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/SRR453110_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453111

echo "* Start fastqc [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/SRR453111_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/SRR453111_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453112

echo "* Start fastqc [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/SRR453112_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/SRR453112_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453113

echo "* Start fastqc [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/SRR453113_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/SRR453113_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453114

echo "* Start fastqc [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/SRR453114_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/SRR453114_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453115

echo "* Start fastqc [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/SRR453115_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/SRR453115_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453109

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed

echo "* Start scythe [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/SRR453109_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/SRR453109_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453110

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed

echo "* Start scythe [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/SRR453110_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/SRR453110_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453111

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed

echo "* Start scythe [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/SRR453111_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/SRR453111_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453112

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed

echo "* Start scythe [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/SRR453112_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/SRR453112_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453113

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed

echo "* Start scythe [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/SRR453113_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/SRR453113_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453114

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed

echo "* Start scythe [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/SRR453114_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/SRR453114_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453115

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed

echo "* Start scythe [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/SRR453115_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/SRR453115_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453109
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed

echo "* Start sickle [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453110
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed

echo "* Start sickle [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453111
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed

echo "* Start sickle [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453112
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed

echo "* Start sickle [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453113
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed

echo "* Start sickle [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453114
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed

echo "* Start sickle [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453115
cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed

echo "* Start sickle [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453109

echo "* Start fastqc [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453110

echo "* Start fastqc [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453111

echo "* Start fastqc [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453112

echo "* Start fastqc [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453113

echo "* Start fastqc [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453114

echo "* Start fastqc [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453115

echo "* Start fastqc [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


