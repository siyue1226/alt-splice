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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070411
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411

echo "* Start srr_dump [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070411.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070428
mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428

echo "* Start srr_dump [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070428.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070411

echo "* Start fastqc [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/SRR070411_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/SRR070411_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070428

echo "* Start fastqc [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/SRR070428_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/SRR070428_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070411

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed

echo "* Start scythe [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/SRR070411_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/SRR070411_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070428

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed

echo "* Start scythe [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/SRR070428_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/SRR070428_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070411
cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed

echo "* Start sickle [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070428
cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed

echo "* Start sickle [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070411

echo "* Start fastqc [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070428

echo "* Start fastqc [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


