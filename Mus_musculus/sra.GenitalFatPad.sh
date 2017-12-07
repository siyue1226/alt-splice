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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453126
mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126

echo "* Start srr_dump [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453126.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453127
mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127

echo "* Start srr_dump [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453127.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453128
mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128

echo "* Start srr_dump [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453128.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453129
mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129

echo "* Start srr_dump [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453129.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453126

echo "* Start fastqc [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/SRR453126_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/SRR453126_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453127

echo "* Start fastqc [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/SRR453127_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/SRR453127_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453128

echo "* Start fastqc [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/SRR453128_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/SRR453128_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453129

echo "* Start fastqc [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/SRR453129_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/SRR453129_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453126

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed

echo "* Start scythe [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/SRR453126_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/SRR453126_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453127

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed

echo "* Start scythe [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/SRR453127_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/SRR453127_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453128

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed

echo "* Start scythe [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/SRR453128_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/SRR453128_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453129

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed

echo "* Start scythe [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/SRR453129_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/SRR453129_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453126
cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed

echo "* Start sickle [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453127
cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed

echo "* Start sickle [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453128
cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed

echo "* Start sickle [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453129
cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed

echo "* Start sickle [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453126

echo "* Start fastqc [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453127

echo "* Start fastqc [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453128

echo "* Start fastqc [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453129

echo "* Start fastqc [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


