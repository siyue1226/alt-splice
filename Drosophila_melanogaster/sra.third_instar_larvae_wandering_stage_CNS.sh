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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070409
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409

echo "* Start srr_dump [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070409.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070409 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070410
mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410

echo "* Start srr_dump [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070410.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070410 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070409

echo "* Start fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/SRR070409_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/SRR070409_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070409 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070410

echo "* Start fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/SRR070410_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/SRR070410_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070410 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070409

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/SRR070409_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/SRR070409_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070409 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070410

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed

echo "* Start scythe [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/SRR070410_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/SRR070410_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070410 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070409
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070409 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070410
cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed

echo "* Start sickle [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070410 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070409

echo "* Start fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/trimmed/SRR070409_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070409 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070409] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070410

echo "* Start fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/trimmed/SRR070410_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_CNS SRR070410 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [third_instar_larvae_wandering_stage_CNS] [SRR070410] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


