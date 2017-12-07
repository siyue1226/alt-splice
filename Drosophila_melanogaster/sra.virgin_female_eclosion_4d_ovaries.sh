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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070396
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396

echo "* Start srr_dump [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070396.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070417
mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417

echo "* Start srr_dump [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070417.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070396

echo "* Start fastqc [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/SRR070396_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/SRR070396_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070417

echo "* Start fastqc [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/SRR070417_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/SRR070417_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070396

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed

echo "* Start scythe [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/SRR070396_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/SRR070396_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070417

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed

echo "* Start scythe [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/SRR070417_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/SRR070417_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070396
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed

echo "* Start sickle [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070417
cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed

echo "* Start sickle [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070396

echo "* Start fastqc [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070417

echo "* Start fastqc [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


