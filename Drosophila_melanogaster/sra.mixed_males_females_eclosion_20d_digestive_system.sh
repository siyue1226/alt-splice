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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070403
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403

echo "* Start srr_dump [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070403.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR111883
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883

echo "* Start srr_dump [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR111883.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070403

echo "* Start fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/SRR070403_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/SRR070403_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111883

echo "* Start fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/SRR111883_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/SRR111883_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070403

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed

echo "* Start scythe [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/SRR070403_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/SRR070403_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR111883

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed

echo "* Start scythe [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/SRR111883_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/SRR111883_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070403
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed

echo "* Start sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR111883
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed

echo "* Start sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070403

echo "* Start fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR111883

echo "* Start fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


