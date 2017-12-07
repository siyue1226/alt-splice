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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070426
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426

echo "* Start srr_dump [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070426.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR100269
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269

echo "* Start srr_dump [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR100269.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070426

echo "* Start fastqc [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/SRR070426_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/SRR070426_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100269

echo "* Start fastqc [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/SRR100269_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/SRR100269_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070426

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/SRR070426_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/SRR070426_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR100269

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/SRR100269_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/SRR100269_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070426
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR100269
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070426

echo "* Start fastqc [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR100269

echo "* Start fastqc [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


