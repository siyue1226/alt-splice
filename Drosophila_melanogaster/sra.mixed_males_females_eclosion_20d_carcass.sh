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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070404
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404

echo "* Start srr_dump [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070404.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070391
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391

echo "* Start srr_dump [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070391.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070404

echo "* Start fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/SRR070404_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/SRR070404_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070391

echo "* Start fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/SRR070391_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/SRR070391_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070404

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed

echo "* Start scythe [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/SRR070404_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/SRR070404_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070391

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed

echo "* Start scythe [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/SRR070391_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/SRR070391_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070404
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed

echo "* Start sickle [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070391
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed

echo "* Start sickle [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070404

echo "* Start fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070391

echo "* Start fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


