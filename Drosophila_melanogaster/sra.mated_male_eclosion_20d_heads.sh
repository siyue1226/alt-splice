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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070421
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421

echo "* Start srr_dump [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070421.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070424
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424

echo "* Start srr_dump [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070424.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070421

echo "* Start fastqc [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/SRR070421_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/SRR070421_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070424

echo "* Start fastqc [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/SRR070424_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/SRR070424_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070421

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed

echo "* Start scythe [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/SRR070421_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/SRR070421_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070424

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed

echo "* Start scythe [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/SRR070424_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/SRR070424_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070421
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed

echo "* Start sickle [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070424
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed

echo "* Start sickle [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070421

echo "* Start fastqc [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070424

echo "* Start fastqc [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


