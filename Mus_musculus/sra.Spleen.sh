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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Spleen ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453160
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160

echo "* Start srr_dump [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453160.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453160 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453161
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161

echo "* Start srr_dump [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453161.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453161 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453162
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162

echo "* Start srr_dump [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453162.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453162 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453163
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163

echo "* Start srr_dump [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453163.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453163 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453164
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164

echo "* Start srr_dump [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453164.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453164 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453165
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165

echo "* Start srr_dump [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453165.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Spleen SRR453165 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453160

echo "* Start fastqc [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/SRR453160_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/SRR453160_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453160 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453161

echo "* Start fastqc [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/SRR453161_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/SRR453161_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453161 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453162

echo "* Start fastqc [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/SRR453162_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/SRR453162_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453162 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453163

echo "* Start fastqc [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/SRR453163_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/SRR453163_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453163 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453164

echo "* Start fastqc [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/SRR453164_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/SRR453164_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453164 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453165

echo "* Start fastqc [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/SRR453165_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/SRR453165_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453165 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453160

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed

echo "* Start scythe [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/SRR453160_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/SRR453160_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453160 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453161

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed

echo "* Start scythe [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/SRR453161_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/SRR453161_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453161 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453162

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed

echo "* Start scythe [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/SRR453162_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/SRR453162_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453162 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453163

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed

echo "* Start scythe [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/SRR453163_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/SRR453163_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453163 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453164

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed

echo "* Start scythe [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/SRR453164_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/SRR453164_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453164 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453165

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed

echo "* Start scythe [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/SRR453165_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/SRR453165_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453165 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453160
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed

echo "* Start sickle [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453160 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453161
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed

echo "* Start sickle [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453161 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453162
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed

echo "* Start sickle [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453162 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453163
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed

echo "* Start sickle [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453163 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453164
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed

echo "* Start sickle [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453164 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453165
cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed

echo "* Start sickle [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453165 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453160

echo "* Start fastqc [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453160 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453161

echo "* Start fastqc [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453161 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453162

echo "* Start fastqc [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453162 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453163

echo "* Start fastqc [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453163 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453164

echo "* Start fastqc [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453164 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453165

echo "* Start fastqc [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Spleen SRR453165 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


