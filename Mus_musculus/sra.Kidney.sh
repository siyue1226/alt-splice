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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Kidney ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453144
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144

echo "* Start srr_dump [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453144.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453144 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453145
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145

echo "* Start srr_dump [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453145.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453145 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453146
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146

echo "* Start srr_dump [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453146.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453146 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453147
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147

echo "* Start srr_dump [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453147.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453147 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453148
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148

echo "* Start srr_dump [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453148.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453148 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453149
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149

echo "* Start srr_dump [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453149.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Kidney SRR453149 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453144

echo "* Start fastqc [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/SRR453144_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/SRR453144_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453144 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453145

echo "* Start fastqc [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/SRR453145_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/SRR453145_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453145 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453146

echo "* Start fastqc [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/SRR453146_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/SRR453146_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453146 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453147

echo "* Start fastqc [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/SRR453147_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/SRR453147_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453147 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453148

echo "* Start fastqc [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/SRR453148_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/SRR453148_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453148 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453149

echo "* Start fastqc [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/SRR453149_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/SRR453149_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453149 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453144

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed

echo "* Start scythe [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/SRR453144_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/SRR453144_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453144 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453145

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed

echo "* Start scythe [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/SRR453145_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/SRR453145_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453145 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453146

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed

echo "* Start scythe [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/SRR453146_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/SRR453146_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453146 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453147

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed

echo "* Start scythe [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/SRR453147_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/SRR453147_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453147 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453148

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed

echo "* Start scythe [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/SRR453148_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/SRR453148_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453148 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453149

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed

echo "* Start scythe [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/SRR453149_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/SRR453149_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453149 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453144
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed

echo "* Start sickle [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453144 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453145
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed

echo "* Start sickle [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453145 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453146
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed

echo "* Start sickle [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453146 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453147
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed

echo "* Start sickle [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453147 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453148
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed

echo "* Start sickle [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453148 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453149
cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed

echo "* Start sickle [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453149 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453144

echo "* Start fastqc [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453144 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453145

echo "* Start fastqc [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453145 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453146

echo "* Start fastqc [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453146 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453147

echo "* Start fastqc [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453147 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453148

echo "* Start fastqc [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453148 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453149

echo "* Start fastqc [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Kidney SRR453149 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


