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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070407
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407

echo "* Start srr_dump [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070407.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070425
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425

echo "* Start srr_dump [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070425.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070407

echo "* Start fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/SRR070407_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/SRR070407_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070425

echo "* Start fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/SRR070425_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/SRR070425_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070407

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/SRR070407_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/SRR070407_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070425

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/SRR070425_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/SRR070425_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070407
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070425
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070407

echo "* Start fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070425

echo "* Start fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


