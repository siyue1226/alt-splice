#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/rice_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa
# samtools faidx /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa

if [ -d /home/wangq/data/rna-seq/rice_trans/process/ ];
then
    rm -fr /home/wangq/data/rna-seq/rice_trans/process/ ;
fi;
mkdir /home/wangq/data/rna-seq/rice_trans/process/

#----------------------------#
# srr dump
#----------------------------#
# lane SRR352189
mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352189

echo "* Start srr_dump [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/rice_trans/SRR352189.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/rice_trans/process/SRR352189 \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` OSN_AC SRR352189 [fastq dump] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End srr_dump [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/srr_dump.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR352189

echo "* Start fastqc [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/SRR352189_1.fastq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/SRR352189_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AC SRR352189 [fastqc] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End fastqc [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log



#----------------------------#
# scythe
#----------------------------#
# lane SRR352189

if [ ! -d /home/wangq/data/rna-seq/rice_trans/process/SRR352189 ];
then
    mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352189;
fi;

if [ ! -d /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed ;
fi;

cd /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed

echo "* Start scythe [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/SRR352189_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/rice_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/SRR352189_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/rice_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` OSN_AC SRR352189 [scythe] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End scythe [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/scythe.log



#----------------------------#
# sickle
#----------------------------#
# lane SRR352189
cd /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed

echo "* Start sickle [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_1.sickle.fq \
    -p /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_2.sickle.fq \
    -s /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AC SRR352189 [sickle] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End sickle [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log

find /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/sickle.log



#----------------------------#
# fastqc
#----------------------------#
# lane SRR352189

echo "* Start fastqc [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/trimmed/SRR352189_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` OSN_AC SRR352189 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End fastqc [OSN_AC] [SRR352189] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/fastqc.log



