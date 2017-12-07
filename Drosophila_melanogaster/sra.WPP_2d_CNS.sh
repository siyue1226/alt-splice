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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070412
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412

echo "* Start srr_dump [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070412.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100271
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271

echo "* Start srr_dump [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100271.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070412

echo "* Start fastqc [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/SRR070412_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/SRR070412_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100271

echo "* Start fastqc [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/SRR100271_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/SRR100271_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070412

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed

echo "* Start scythe [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/SRR070412_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/SRR070412_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100271

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed

echo "* Start scythe [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/SRR100271_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/SRR100271_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070412
cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed

echo "* Start sickle [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100271
cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed

echo "* Start sickle [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070412

echo "* Start fastqc [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100271

echo "* Start fastqc [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


