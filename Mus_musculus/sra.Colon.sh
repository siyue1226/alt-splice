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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Colon ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Colon ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453166
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166

echo "* Start srr_dump [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453166.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453166 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453167
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167

echo "* Start srr_dump [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453167.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453167 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453168
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168

echo "* Start srr_dump [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453168.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453168 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453169
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169

echo "* Start srr_dump [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453169.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453169 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453170
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170

echo "* Start srr_dump [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453170.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453170 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453171
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171

echo "* Start srr_dump [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453171.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Colon SRR453171 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453166

echo "* Start fastqc [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/SRR453166_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/SRR453166_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453166 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453167

echo "* Start fastqc [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/SRR453167_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/SRR453167_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453167 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453168

echo "* Start fastqc [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/SRR453168_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/SRR453168_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453168 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453169

echo "* Start fastqc [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/SRR453169_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/SRR453169_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453169 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453170

echo "* Start fastqc [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/SRR453170_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/SRR453170_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453170 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453171

echo "* Start fastqc [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/SRR453171_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/SRR453171_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453171 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453166

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed

echo "* Start scythe [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/SRR453166_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/SRR453166_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453166 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453167

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed

echo "* Start scythe [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/SRR453167_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/SRR453167_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453167 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453168

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed

echo "* Start scythe [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/SRR453168_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/SRR453168_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453168 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453169

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed

echo "* Start scythe [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/SRR453169_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/SRR453169_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453169 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453170

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed

echo "* Start scythe [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/SRR453170_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/SRR453170_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453170 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453171

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed

echo "* Start scythe [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/SRR453171_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/SRR453171_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453171 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453166
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed

echo "* Start sickle [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453166 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453167
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed

echo "* Start sickle [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453167 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453168
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed

echo "* Start sickle [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453168 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453169
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed

echo "* Start sickle [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453169 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453170
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed

echo "* Start sickle [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453170 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453171
cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed

echo "* Start sickle [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453171 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453166

echo "* Start fastqc [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453166 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453167

echo "* Start fastqc [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453167 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453168

echo "* Start fastqc [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453168 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453169

echo "* Start fastqc [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453169 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453170

echo "* Start fastqc [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453170 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453171

echo "* Start fastqc [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Colon SRR453171 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


