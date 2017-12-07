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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070400
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400

echo "* Start srr_dump [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070400.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070416
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416

echo "* Start srr_dump [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070416.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070400

echo "* Start fastqc [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/SRR070400_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/SRR070400_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070416

echo "* Start fastqc [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/SRR070416_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/SRR070416_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070400

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed

echo "* Start scythe [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/SRR070400_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/SRR070400_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070416

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed

echo "* Start scythe [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/SRR070416_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/SRR070416_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070400
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed

echo "* Start sickle [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070416
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed

echo "* Start sickle [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070400

echo "* Start fastqc [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070416

echo "* Start fastqc [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


