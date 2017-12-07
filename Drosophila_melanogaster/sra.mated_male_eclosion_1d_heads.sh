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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070432
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432

echo "* Start srr_dump [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070432.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070433
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433

echo "* Start srr_dump [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070433.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100280
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280

echo "* Start srr_dump [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100280.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070432

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/SRR070432_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/SRR070432_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070433

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/SRR070433_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/SRR070433_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100280

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/SRR100280_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/SRR100280_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070432

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed

echo "* Start scythe [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/SRR070432_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/SRR070432_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070433

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed

echo "* Start scythe [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/SRR070433_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/SRR070433_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100280

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed

echo "* Start scythe [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/SRR100280_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/SRR100280_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070432
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed

echo "* Start sickle [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070433
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed

echo "* Start sickle [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100280
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed

echo "* Start sickle [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070432

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070433

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100280

echo "* Start fastqc [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


