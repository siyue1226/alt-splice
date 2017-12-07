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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453130
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130

echo "* Start srr_dump [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453130.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453131
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131

echo "* Start srr_dump [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453131.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453132
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132

echo "* Start srr_dump [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453132.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453133
mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133

echo "* Start srr_dump [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453133.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453130

echo "* Start fastqc [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/SRR453130_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/SRR453130_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453131

echo "* Start fastqc [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/SRR453131_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/SRR453131_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453132

echo "* Start fastqc [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/SRR453132_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/SRR453132_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453133

echo "* Start fastqc [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/SRR453133_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/SRR453133_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453130

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed

echo "* Start scythe [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/SRR453130_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/SRR453130_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453131

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed

echo "* Start scythe [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/SRR453131_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/SRR453131_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453132

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed

echo "* Start scythe [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/SRR453132_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/SRR453132_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453133

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed

echo "* Start scythe [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/SRR453133_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/SRR453133_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453130
cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed

echo "* Start sickle [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453131
cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed

echo "* Start sickle [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453132
cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed

echo "* Start sickle [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453133
cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed

echo "* Start sickle [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453130

echo "* Start fastqc [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453131

echo "* Start fastqc [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453132

echo "* Start fastqc [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453133

echo "* Start fastqc [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


