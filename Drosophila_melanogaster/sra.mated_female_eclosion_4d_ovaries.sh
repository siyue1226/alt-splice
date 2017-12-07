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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070431
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431

echo "* Start srr_dump [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070431.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100277
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277

echo "* Start srr_dump [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100277.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100283
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283

echo "* Start srr_dump [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100283.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070431

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/SRR070431_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/SRR070431_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100277

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/SRR100277_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/SRR100277_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100283

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/SRR100283_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/SRR100283_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070431

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed

echo "* Start scythe [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/SRR070431_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/SRR070431_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100277

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed

echo "* Start scythe [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/SRR100277_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/SRR100277_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100283

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed

echo "* Start scythe [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/SRR100283_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/SRR100283_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070431
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed

echo "* Start sickle [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100277
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed

echo "* Start sickle [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100283
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed

echo "* Start sickle [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070431

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100277

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100283

echo "* Start fastqc [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


