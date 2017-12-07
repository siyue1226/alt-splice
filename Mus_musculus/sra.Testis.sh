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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Testis ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Testis ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453140
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140

echo "* Start srr_dump [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453140.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Testis SRR453140 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453141
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141

echo "* Start srr_dump [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453141.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Testis SRR453141 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453142
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142

echo "* Start srr_dump [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453142.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Testis SRR453142 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453143
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143

echo "* Start srr_dump [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453143.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Testis SRR453143 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453140

echo "* Start fastqc [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/SRR453140_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/SRR453140_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453140 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453141

echo "* Start fastqc [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/SRR453141_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/SRR453141_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453141 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453142

echo "* Start fastqc [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/SRR453142_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/SRR453142_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453142 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453143

echo "* Start fastqc [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/SRR453143_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/SRR453143_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453143 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453140

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed

echo "* Start scythe [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/SRR453140_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/SRR453140_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453140 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453141

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed

echo "* Start scythe [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/SRR453141_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/SRR453141_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453141 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453142

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed

echo "* Start scythe [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/SRR453142_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/SRR453142_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453142 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453143

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed

echo "* Start scythe [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/SRR453143_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/SRR453143_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453143 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453140
cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed

echo "* Start sickle [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453140 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453141
cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed

echo "* Start sickle [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453141 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453142
cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed

echo "* Start sickle [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453142 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453143
cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed

echo "* Start sickle [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453143 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453140

echo "* Start fastqc [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453140 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453141

echo "* Start fastqc [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453141 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453142

echo "* Start fastqc [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453142 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453143

echo "* Start fastqc [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Testis SRR453143 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


