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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Ovary ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453077
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077

echo "* Start srr_dump [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453077.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453077 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453078
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078

echo "* Start srr_dump [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453078.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453078 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453079
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079

echo "* Start srr_dump [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453079.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453079 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453080
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080

echo "* Start srr_dump [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453080.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453080 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453081
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081

echo "* Start srr_dump [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453081.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453081 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453082
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082

echo "* Start srr_dump [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453082.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453082 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453083
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083

echo "* Start srr_dump [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453083.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453083 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453084
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084

echo "* Start srr_dump [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453084.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453084 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453085
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085

echo "* Start srr_dump [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453085.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453085 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453086
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086

echo "* Start srr_dump [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453086.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Ovary SRR453086 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453077

echo "* Start fastqc [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/SRR453077_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/SRR453077_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453077 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453078

echo "* Start fastqc [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/SRR453078_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/SRR453078_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453078 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453079

echo "* Start fastqc [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/SRR453079_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/SRR453079_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453079 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453080

echo "* Start fastqc [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/SRR453080_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/SRR453080_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453080 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453081

echo "* Start fastqc [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/SRR453081_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/SRR453081_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453081 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453082

echo "* Start fastqc [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/SRR453082_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/SRR453082_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453082 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453083

echo "* Start fastqc [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/SRR453083_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/SRR453083_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453083 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453084

echo "* Start fastqc [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/SRR453084_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/SRR453084_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453084 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453085

echo "* Start fastqc [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/SRR453085_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/SRR453085_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453085 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453086

echo "* Start fastqc [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/SRR453086_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/SRR453086_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453086 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453077

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed

echo "* Start scythe [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/SRR453077_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/SRR453077_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453077 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453078

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed

echo "* Start scythe [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/SRR453078_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/SRR453078_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453078 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453079

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed

echo "* Start scythe [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/SRR453079_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/SRR453079_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453079 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453080

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed

echo "* Start scythe [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/SRR453080_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/SRR453080_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453080 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453081

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed

echo "* Start scythe [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/SRR453081_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/SRR453081_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453081 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453082

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed

echo "* Start scythe [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/SRR453082_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/SRR453082_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453082 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453083

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed

echo "* Start scythe [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/SRR453083_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/SRR453083_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453083 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453084

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed

echo "* Start scythe [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/SRR453084_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/SRR453084_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453084 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453085

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed

echo "* Start scythe [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/SRR453085_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/SRR453085_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453085 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453086

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed

echo "* Start scythe [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/SRR453086_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/SRR453086_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453086 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453077
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed

echo "* Start sickle [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453077 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453078
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed

echo "* Start sickle [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453078 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453079
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed

echo "* Start sickle [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453079 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453080
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed

echo "* Start sickle [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453080 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453081
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed

echo "* Start sickle [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453081 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453082
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed

echo "* Start sickle [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453082 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453083
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed

echo "* Start sickle [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453083 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453084
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed

echo "* Start sickle [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453084 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453085
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed

echo "* Start sickle [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453085 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453086
cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed

echo "* Start sickle [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453086 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453077

echo "* Start fastqc [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453077 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453078

echo "* Start fastqc [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453078 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453079

echo "* Start fastqc [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453079 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453080

echo "* Start fastqc [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453080 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453081

echo "* Start fastqc [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453081 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453082

echo "* Start fastqc [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453082 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453083

echo "* Start fastqc [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453083 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453084

echo "* Start fastqc [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453084 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453085

echo "* Start fastqc [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453085 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453086

echo "* Start fastqc [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Ovary SRR453086 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


