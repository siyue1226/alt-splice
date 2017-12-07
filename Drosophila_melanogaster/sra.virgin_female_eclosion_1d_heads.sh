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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070436
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436

echo "* Start srr_dump [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070436.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070437
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437

echo "* Start srr_dump [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070437.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100281
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281

echo "* Start srr_dump [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100281.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070436

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/SRR070436_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/SRR070436_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070437

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/SRR070437_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/SRR070437_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100281

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/SRR100281_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/SRR100281_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070436

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed

echo "* Start scythe [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/SRR070436_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/SRR070436_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070437

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed

echo "* Start scythe [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/SRR070437_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/SRR070437_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100281

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed

echo "* Start scythe [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/SRR100281_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/SRR100281_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070436
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed

echo "* Start sickle [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070437
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed

echo "* Start sickle [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100281
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed

echo "* Start sickle [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070436

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070437

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100281

echo "* Start fastqc [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


