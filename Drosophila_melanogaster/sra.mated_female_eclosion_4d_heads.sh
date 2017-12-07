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

if [ -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads ];
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads ;
fi;
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads

#----------------------------#
# srr dump
#----------------------------#
# lane SRR070414
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414

echo "* Start srr_dump [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070414.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# lane SRR070415
mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415

echo "* Start srr_dump [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/dmel_trans/SRP001696/SRR070415.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415 \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [fastq dump] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End srr_dump [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070414

echo "* Start fastqc [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/SRR070414_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/SRR070414_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070415

echo "* Start fastqc [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/SRR070415_1.fastq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/SRR070415_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [fastqc] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR070414

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed

echo "* Start scythe [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/SRR070414_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/SRR070414_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# lane SRR070415

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415 ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415;
fi;

if [ ! -d /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed ;
fi;

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed

echo "* Start scythe [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/SRR070415_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/SRR070415_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/dmel_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [scythe] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End scythe [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR070414
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed

echo "* Start sickle [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# lane SRR070415
cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed

echo "* Start sickle [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.sickle.fq \
    -p /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.sickle.fq \
    -s /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End sickle [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log

find /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR070414

echo "* Start fastqc [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# lane SRR070415

echo "* Start fastqc [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End fastqc [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/fastqc.log


