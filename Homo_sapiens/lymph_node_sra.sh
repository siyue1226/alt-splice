#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq1/data/rna2dna/bodymap2

### index reference genome
# bwa index -a bwtsw /home/wangq1/data/rna2dna/bodymap2/ref/human.37.fa
# samtools faidx /home/wangq1/data/rna2dna/bodymap2/ref/human.37.fa

if [ -d /home/wangq1/data/rna2dna/bodymap2/process/lymph_node ];
then
    rm -fr /home/wangq1/data/rna2dna/bodymap2/process/lymph_node ;
fi;
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node

#----------------------------#
# srr dump
#----------------------------#
# lane ERR030897
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897

echo "* Start srr_dump [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (single end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030897.sra \
    --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# lane ERR030878
mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878

echo "* Start srr_dump [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log

# sra to fastq (pair end)
fastq-dump /home/wangq1/data/rna2dna/bodymap2/ERP000546/ERR030878.sra \
    --split-files --gzip -O /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878 \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [fastq dump] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End srr_dump [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030897

echo "* Start fastqc [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/ERR030897.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030878

echo "* Start fastqc [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/ERR030878_1.fastq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/ERR030878_2.fastq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [fastqc] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane ERR030897

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed

echo "* Start scythe [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (single end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/ERR030897.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.scythe.fq.gz


[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# lane ERR030878

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878 ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878;
fi;

if [ ! -d /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed  ];
then
    mkdir /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed ;
fi;

cd /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed

echo "* Start scythe [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log

# scythe (pair end)
scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/ERR030878_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.scythe.fq.gz

scythe \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/ERR030878_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq1/data/rna2dna/bodymap2/ref/illumina_adapters.fa \
    -m /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [scythe] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End scythe [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane ERR030897
cd /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed

echo "* Start sickle [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (single end)
sickle se \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.sickle.fq
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# lane ERR030878
cd /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed

echo "* Start sickle [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

# sickle (pair end)
sickle pe \
    -t sanger \
    -q 20 \
    -l 20 \
    -f /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.scythe.fq.gz \
    -r /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.scythe.fq.gz \
    -o /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.sickle.fq \
    -p /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.sickle.fq \
    -s /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_single.sickle.fq \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End sickle [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log

find /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane ERR030897

echo "* Start fastqc [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (single end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lymph_node] [ERR030897] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# lane ERR030878

echo "* Start fastqc [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.sickle.fq.gz \
    /home/wangq1/data/rna2dna/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [fastqc.sickle] failed >> /home/wangq1/data/rna2dna/bodymap2/fail.log && exit 255
echo "* End fastqc [lymph_node] [ERR030878] `date`" | tee -a /home/wangq1/data/rna2dna/bodymap2/log/fastqc.log


