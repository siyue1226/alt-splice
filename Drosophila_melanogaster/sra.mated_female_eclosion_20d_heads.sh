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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070420
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420

echo "* Start srr_dump [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070420.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100274
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274

echo "* Start srr_dump [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100274.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR116383
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383

echo "* Start srr_dump [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR116383.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111882
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882

echo "* Start srr_dump [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111882.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070420

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/SRR070420_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/SRR070420_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100274

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/SRR100274_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/SRR100274_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR116383

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/SRR116383_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/SRR116383_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111882

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/SRR111882_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/SRR111882_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070420

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed

echo "* Start scythe [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/SRR070420_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/SRR070420_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100274

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed

echo "* Start scythe [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/SRR100274_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/SRR100274_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR116383

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed

echo "* Start scythe [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/SRR116383_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/SRR116383_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111882

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed

echo "* Start scythe [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/SRR111882_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/SRR111882_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070420
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed

echo "* Start sickle [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100274
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed

echo "* Start sickle [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR116383
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed

echo "* Start sickle [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111882
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed

echo "* Start sickle [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070420

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100274

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR116383

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111882

echo "* Start fastqc [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


