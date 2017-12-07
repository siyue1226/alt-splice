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

if [ -d /home/wangq/data/rna-seq/mouse_trans/process/Heart ];
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/Heart ;
fi;
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart

#----------------------------#
# srr dump
#----------------------------#
# lane SRR453172
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172

echo "* Start srr_dump [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453172.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Heart SRR453172 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453173
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173

echo "* Start srr_dump [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453173.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Heart SRR453173 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453174
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174

echo "* Start srr_dump [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453174.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Heart SRR453174 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# lane SRR453175
mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175

echo "* Start srr_dump [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log

# sra to fastq (pair end)
/home/wangq/share/sratoolkit/fastq-dump /home/wangq/data/rna-seq/mouse_trans/SRP012040/SRR453175.sra \
    --split-files --gzip -O /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175 \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log ; ( exit ${PIPESTATUS} )

[ $? -ne 0 ] && echo `date` Heart SRR453175 [fastq dump] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End srr_dump [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/srr_dump.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453172

echo "* Start fastqc [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/SRR453172_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/SRR453172_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453172 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453173

echo "* Start fastqc [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/SRR453173_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/SRR453173_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453173 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453174

echo "* Start fastqc [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/SRR453174_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/SRR453174_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453174 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453175

echo "* Start fastqc [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/SRR453175_1.fastq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/SRR453175_2.fastq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453175 [fastqc] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


#----------------------------#
# scythe
#----------------------------#
# lane SRR453172

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed

echo "* Start scythe [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/SRR453172_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/SRR453172_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453172 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453173

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed

echo "* Start scythe [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/SRR453173_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/SRR453173_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453173 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453174

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed

echo "* Start scythe [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/SRR453174_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/SRR453174_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453174 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453175

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed

echo "* Start scythe [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/SRR453175_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/SRR453175_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453175 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453172
cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed

echo "* Start sickle [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453172 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453173
cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed

echo "* Start sickle [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453173 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453174
cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed

echo "* Start sickle [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453174 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453175
cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed

echo "* Start sickle [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453175 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453172

echo "* Start fastqc [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453172 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453173

echo "* Start fastqc [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453173 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453174

echo "* Start fastqc [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453174 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453175

echo "* Start fastqc [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Heart SRR453175 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


