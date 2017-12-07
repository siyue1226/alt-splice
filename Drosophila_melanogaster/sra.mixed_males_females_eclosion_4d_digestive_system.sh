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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070401
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401

echo "* Start srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070401.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111878
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878

echo "* Start srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111878.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111879
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879

echo "* Start srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111879.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070401

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/SRR070401_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/SRR070401_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111878

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/SRR111878_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/SRR111878_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111879

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/SRR111879_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/SRR111879_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070401

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed

echo "* Start scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/SRR070401_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/SRR070401_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111878

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed

echo "* Start scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/SRR111878_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/SRR111878_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111879

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed

echo "* Start scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/SRR111879_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/SRR111879_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070401
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed

echo "* Start sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111878
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed

echo "* Start sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111879
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed

echo "* Start sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070401

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111878

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111879

echo "* Start fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


