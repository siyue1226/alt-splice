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

#----------------------------#
# scythe
#----------------------------#
# lane SRR453150

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed

echo "* Start scythe [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/SRR453150_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/SRR453150_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453150 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453151

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed

echo "* Start scythe [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/SRR453151_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/SRR453151_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453151 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453152

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed

echo "* Start scythe [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/SRR453152_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/SRR453152_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453152 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453153

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed

echo "* Start scythe [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/SRR453153_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/SRR453153_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453153 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453154

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed

echo "* Start scythe [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/SRR453154_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/SRR453154_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453154 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# lane SRR453155

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155 ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155;
fi;

if [ ! -d /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed  ];
then
    mkdir /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed ;
fi;

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed

echo "* Start scythe [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log

# scythe (pair end)
/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/SRR453155_1.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.scythe.fq.gz

/home/wangq/bin/scythe \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/SRR453155_2.fastq.gz \
    -q sanger \
    -M 20 \
    -a /home/wangq/data/rna-seq/mouse_trans/ref/illumina_adapters.fa \
    -m /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.matches.txt \
    --quiet \
    | gzip -c --fast > /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.scythe.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453155 [scythe] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End scythe [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/scythe.log


#----------------------------#
# sickle
#----------------------------#
# lane SRR453150
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed

echo "* Start sickle [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453150 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453151
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed

echo "* Start sickle [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453151 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453152
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed

echo "* Start sickle [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453152 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453153
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed

echo "* Start sickle [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453153 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453154
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed

echo "* Start sickle [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453154 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# lane SRR453155
cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed

echo "* Start sickle [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

# sickle (pair end)
/home/wangq/bin/sickle pe \
    -t sanger -l 20 -q 20 \
    -f /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.scythe.fq.gz \
    -r /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.scythe.fq.gz \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.sickle.fq \
    -p /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.sickle.fq \
    -s /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_single.sickle.fq \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453155 [sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End sickle [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log

find /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/ -type f -name "*.sickle.fq" | parallel -j 8 gzip --best
echo "* Gzip sickle [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/sickle.log


#----------------------------#
# fastqc
#----------------------------#
# lane SRR453150

echo "* Start fastqc [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453150 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453151

echo "* Start fastqc [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453151 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453152

echo "* Start fastqc [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453152 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453153

echo "* Start fastqc [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453153 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453154

echo "* Start fastqc [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453154 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# lane SRR453155

echo "* Start fastqc [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log

# fastqc (pair end)
/home/wangq/bin/fastqc -t 8 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_single.sickle.fq.gz \
    2>&1 | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log ; ( exit ${PIPESTATUS} )


[ $? -ne 0 ] && echo `date` Liver SRR453155 [fastqc.sickle] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End fastqc [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/fastqc.log


