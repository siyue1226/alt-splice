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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070434
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434

echo "* Start srr_dump [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070434.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070435
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435

echo "* Start srr_dump [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070435.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100279
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279

echo "* Start srr_dump [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100279.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070434

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/SRR070434_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/SRR070434_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070435

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/SRR070435_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/SRR070435_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100279

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/SRR100279_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/SRR100279_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070434

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed

echo "* Start scythe [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/SRR070434_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/SRR070434_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070435

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed

echo "* Start scythe [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/SRR070435_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/SRR070435_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100279

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed

echo "* Start scythe [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/SRR100279_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/SRR100279_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070434
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed

echo "* Start sickle [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070435
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed

echo "* Start sickle [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100279
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed

echo "* Start sickle [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070434

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/trimmed/SRR070434_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070435

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/trimmed/SRR070435_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100279

echo "* Start fastqc [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/trimmed/SRR100279_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


