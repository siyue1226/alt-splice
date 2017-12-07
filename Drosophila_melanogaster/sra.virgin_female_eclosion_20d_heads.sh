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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070388
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388

echo "* Start srr_dump [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070388.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070419
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419

echo "* Start srr_dump [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070419.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100275
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275

echo "* Start srr_dump [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100275.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070388

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/SRR070388_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/SRR070388_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070419

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/SRR070419_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/SRR070419_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100275

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/SRR100275_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/SRR100275_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070388

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed

echo "* Start scythe [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/SRR070388_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/SRR070388_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070419

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed

echo "* Start scythe [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/SRR070419_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/SRR070419_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100275

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed

echo "* Start scythe [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/SRR100275_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/SRR100275_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070388
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed

echo "* Start sickle [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070419
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed

echo "* Start sickle [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100275
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed

echo "* Start sickle [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070388

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070419

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100275

echo "* Start fastqc [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


