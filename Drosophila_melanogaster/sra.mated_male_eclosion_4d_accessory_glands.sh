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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070397
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397

echo "* Start srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070397.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR182358
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358

echo "* Start srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR182358.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR350959
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959

echo "* Start srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR350959.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070397

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/SRR070397_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/SRR070397_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR182358

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/SRR182358_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/SRR182358_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350959

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/SRR350959_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/SRR350959_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070397

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed

echo "* Start scythe [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/SRR070397_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/SRR070397_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR182358

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed

echo "* Start scythe [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/SRR182358_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/SRR182358_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR350959

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed

echo "* Start scythe [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/SRR350959_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/SRR350959_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070397
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed

echo "* Start sickle [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR182358
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed

echo "* Start sickle [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR350959
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed

echo "* Start sickle [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070397

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR182358

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR350959

echo "* Start fastqc [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


