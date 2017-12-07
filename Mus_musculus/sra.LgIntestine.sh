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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453122
mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122

echo "* Start srr_dump [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453122.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453123
mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123

echo "* Start srr_dump [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453123.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453124
mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124

echo "* Start srr_dump [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453124.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453125
mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125

echo "* Start srr_dump [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453125.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453122

echo "* Start fastqc [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/SRR453122_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/SRR453122_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453123

echo "* Start fastqc [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/SRR453123_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/SRR453123_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453124

echo "* Start fastqc [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/SRR453124_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/SRR453124_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453125

echo "* Start fastqc [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/SRR453125_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/SRR453125_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453122

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed

echo "* Start scythe [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/SRR453122_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/SRR453122_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453123

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed

echo "* Start scythe [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/SRR453123_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/SRR453123_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453124

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed

echo "* Start scythe [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/SRR453124_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/SRR453124_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453125

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed

echo "* Start scythe [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/SRR453125_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/SRR453125_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453122
cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed

echo "* Start sickle [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453123
cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed

echo "* Start sickle [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453124
cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed

echo "* Start sickle [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453125
cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed

echo "* Start sickle [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453122

echo "* Start fastqc [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453123

echo "* Start fastqc [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453124

echo "* Start fastqc [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453125

echo "* Start fastqc [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


