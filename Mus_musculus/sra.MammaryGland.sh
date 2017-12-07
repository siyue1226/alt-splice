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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453087
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087

echo "* Start srr_dump [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453087.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453088
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088

echo "* Start srr_dump [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453088.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453089
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089

echo "* Start srr_dump [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453089.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453090
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090

echo "* Start srr_dump [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453090.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453091
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091

echo "* Start srr_dump [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453091.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453092
mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092

echo "* Start srr_dump [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453092.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453087

echo "* Start fastqc [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/SRR453087_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/SRR453087_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453088

echo "* Start fastqc [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/SRR453088_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/SRR453088_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453089

echo "* Start fastqc [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/SRR453089_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/SRR453089_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453090

echo "* Start fastqc [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/SRR453090_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/SRR453090_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453091

echo "* Start fastqc [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/SRR453091_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/SRR453091_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453092

echo "* Start fastqc [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/SRR453092_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/SRR453092_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453087

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed

echo "* Start scythe [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/SRR453087_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/SRR453087_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453088

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed

echo "* Start scythe [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/SRR453088_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/SRR453088_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453089

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed

echo "* Start scythe [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/SRR453089_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/SRR453089_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453090

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed

echo "* Start scythe [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/SRR453090_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/SRR453090_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453091

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed

echo "* Start scythe [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/SRR453091_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/SRR453091_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453092

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed

echo "* Start scythe [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/SRR453092_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/SRR453092_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453087
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed

echo "* Start sickle [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453088
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed

echo "* Start sickle [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453089
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed

echo "* Start sickle [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453090
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed

echo "* Start sickle [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453091
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed

echo "* Start sickle [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453092
cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed

echo "* Start sickle [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453087

echo "* Start fastqc [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453088

echo "* Start fastqc [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453089

echo "* Start fastqc [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453090

echo "* Start fastqc [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453091

echo "* Start fastqc [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453092

echo "* Start fastqc [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


