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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070399
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399

echo "* Start srr_dump [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070399.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070395
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395

echo "* Start srr_dump [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070395.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070399

echo "* Start fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/SRR070399_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/SRR070399_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070395

echo "* Start fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/SRR070395_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/SRR070395_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070399

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed

echo "* Start scythe [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/SRR070399_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/SRR070399_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070395

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed

echo "* Start scythe [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/SRR070395_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/SRR070395_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070399
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed

echo "* Start sickle [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070395
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed

echo "* Start sickle [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070399

echo "* Start fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070395

echo "* Start fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


