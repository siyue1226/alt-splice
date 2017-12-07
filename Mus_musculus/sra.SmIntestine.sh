#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/mouse_trans

### index reference genome
# bwa index -a bwtsw /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa
# samtools faidx /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453099
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099

echo "* Start srr_dump [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453099.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453100
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100

echo "* Start srr_dump [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453100.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453101
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101

echo "* Start srr_dump [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453101.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453102
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102

echo "* Start srr_dump [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453102.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453103
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103

echo "* Start srr_dump [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453103.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453104
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104

echo "* Start srr_dump [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453104.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453105
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105

echo "* Start srr_dump [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453105.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453106
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106

echo "* Start srr_dump [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453106.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453107
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107

echo "* Start srr_dump [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453107.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453108
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108

echo "* Start srr_dump [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453108.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453099

echo "* Start fastqc [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/SRR453099_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/SRR453099_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453100

echo "* Start fastqc [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/SRR453100_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/SRR453100_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453101

echo "* Start fastqc [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/SRR453101_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/SRR453101_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453102

echo "* Start fastqc [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/SRR453102_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/SRR453102_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453103

echo "* Start fastqc [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/SRR453103_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/SRR453103_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453104

echo "* Start fastqc [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/SRR453104_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/SRR453104_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453105

echo "* Start fastqc [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/SRR453105_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/SRR453105_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453106

echo "* Start fastqc [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/SRR453106_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/SRR453106_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453107

echo "* Start fastqc [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/SRR453107_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/SRR453107_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453108

echo "* Start fastqc [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/SRR453108_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/SRR453108_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453099

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed

echo "* Start scythe [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/SRR453099_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/SRR453099_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453100

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed

echo "* Start scythe [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/SRR453100_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/SRR453100_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453101

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed

echo "* Start scythe [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/SRR453101_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/SRR453101_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453102

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed

echo "* Start scythe [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/SRR453102_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/SRR453102_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453103

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed

echo "* Start scythe [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/SRR453103_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/SRR453103_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453104

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed

echo "* Start scythe [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/SRR453104_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/SRR453104_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453105

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed

echo "* Start scythe [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/SRR453105_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/SRR453105_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453106

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed

echo "* Start scythe [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/SRR453106_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/SRR453106_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453107

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed

echo "* Start scythe [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/SRR453107_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/SRR453107_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453108

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed

echo "* Start scythe [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/SRR453108_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/SRR453108_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453099
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed

echo "* Start sickle [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453100
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed

echo "* Start sickle [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453101
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed

echo "* Start sickle [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453102
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed

echo "* Start sickle [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453103
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed

echo "* Start sickle [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453104
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed

echo "* Start sickle [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453105
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed

echo "* Start sickle [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453106
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed

echo "* Start sickle [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453107
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed

echo "* Start sickle [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453108
cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed

echo "* Start sickle [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453099

echo "* Start fastqc [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453100

echo "* Start fastqc [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453101

echo "* Start fastqc [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453102

echo "* Start fastqc [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453103

echo "* Start fastqc [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453104

echo "* Start fastqc [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453105

echo "* Start fastqc [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453106

echo "* Start fastqc [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453107

echo "* Start fastqc [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453108

echo "* Start fastqc [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


