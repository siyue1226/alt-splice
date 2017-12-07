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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070430
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430

echo "* Start srr_dump [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070430.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100278
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278

echo "* Start srr_dump [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100278.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100282
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282

echo "* Start srr_dump [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100282.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070430

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/SRR070430_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/SRR070430_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100278

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/SRR100278_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/SRR100278_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100282

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/SRR100282_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/SRR100282_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070430

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed

echo "* Start scythe [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/SRR070430_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/SRR070430_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100278

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed

echo "* Start scythe [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/SRR100278_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/SRR100278_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100282

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed

echo "* Start scythe [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/SRR100282_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/SRR100282_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070430
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed

echo "* Start sickle [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100278
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed

echo "* Start sickle [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100282
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed

echo "* Start sickle [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070430

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100278

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100282

echo "* Start fastqc [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


