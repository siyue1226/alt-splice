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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Lung ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Lung ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453156
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156

echo "* Start srr_dump [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453156.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Lung SRR453156 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453157
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157

echo "* Start srr_dump [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453157.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Lung SRR453157 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453158
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158

echo "* Start srr_dump [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453158.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Lung SRR453158 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453159
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159

echo "* Start srr_dump [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453159.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Lung SRR453159 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453156

echo "* Start fastqc [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/SRR453156_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/SRR453156_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453156 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453157

echo "* Start fastqc [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/SRR453157_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/SRR453157_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453157 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453158

echo "* Start fastqc [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/SRR453158_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/SRR453158_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453158 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453159

echo "* Start fastqc [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/SRR453159_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/SRR453159_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453159 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453156

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed

echo "* Start scythe [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/SRR453156_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/SRR453156_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453156 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453157

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed

echo "* Start scythe [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/SRR453157_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/SRR453157_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453157 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453158

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed

echo "* Start scythe [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/SRR453158_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/SRR453158_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453158 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453159

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed

echo "* Start scythe [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/SRR453159_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/SRR453159_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453159 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453156
cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed

echo "* Start sickle [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453156 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453157
cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed

echo "* Start sickle [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453157 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453158
cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed

echo "* Start sickle [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453158 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453159
cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed

echo "* Start sickle [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453159 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453156

echo "* Start fastqc [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453156 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453157

echo "* Start fastqc [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453157 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453158

echo "* Start fastqc [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453158 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453159

echo "* Start fastqc [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Lung SRR453159 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


