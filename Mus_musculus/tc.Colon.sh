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
# lane SRR453166

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/

echo "* Start tophat [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/trimmed/SRR453166_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453166 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453167

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/

echo "* Start tophat [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/trimmed/SRR453167_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453167 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453168

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/

echo "* Start tophat [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/trimmed/SRR453168_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453168 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453169

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/

echo "* Start tophat [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/trimmed/SRR453169_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453169 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453170

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/

echo "* Start tophat [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/trimmed/SRR453170_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453170 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453171

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/

echo "* Start tophat [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/trimmed/SRR453171_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Colon SRR453171 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453166

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/

echo "* Start cufflinks [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453166 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453167

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/

echo "* Start cufflinks [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453167 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453168

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/

echo "* Start cufflinks [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453168 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453169

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/

echo "* Start cufflinks [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453169 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453170

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/

echo "* Start cufflinks [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453170 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453171

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/

echo "* Start cufflinks [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453171 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


