#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/mouse_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65

#----------------------------#
# tophat
#----------------------------#
# lane SRR453126

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/

echo "* Start tophat [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/trimmed/SRR453126_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453127

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/

echo "* Start tophat [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/trimmed/SRR453127_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453128

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/

echo "* Start tophat [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/trimmed/SRR453128_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453129

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/

echo "* Start tophat [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/trimmed/SRR453129_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453126

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/

echo "* Start cufflinks [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453127

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/

echo "* Start cufflinks [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453128

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/

echo "* Start cufflinks [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453129

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/

echo "* Start cufflinks [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


