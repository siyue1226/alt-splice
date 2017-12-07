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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070398
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398

echo "* Start srr_dump [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070398.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070398 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070394
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394

echo "* Start srr_dump [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070394.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070394 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070398

echo "* Start fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/SRR070398_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/SRR070398_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070398 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070394

echo "* Start fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/SRR070394_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/SRR070394_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070394 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070398

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed

echo "* Start scythe [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/SRR070398_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/SRR070398_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070398 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070394

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed

echo "* Start scythe [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/SRR070394_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/SRR070394_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070394 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070398
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed

echo "* Start sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070398 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070394
cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed

echo "* Start sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070394 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070398

echo "* Start fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/trimmed/SRR070398_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070398 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070398] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070394

echo "* Start fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/trimmed/SRR070394_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_digestive_system SRR070394 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mixed_males_females_eclosion_1d_digestive_system] [SRR070394] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


