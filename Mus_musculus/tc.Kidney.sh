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
# lane SRR453144

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/

echo "* Start tophat [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/trimmed/SRR453144_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453144 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453145

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/

echo "* Start tophat [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/trimmed/SRR453145_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453145 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453146

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/

echo "* Start tophat [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/trimmed/SRR453146_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453146 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453147

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/

echo "* Start tophat [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/trimmed/SRR453147_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453147 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453148

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/

echo "* Start tophat [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/trimmed/SRR453148_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453148 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453149

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/

echo "* Start tophat [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/trimmed/SRR453149_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Kidney SRR453149 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453144

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/

echo "* Start cufflinks [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453144 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453145

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/

echo "* Start cufflinks [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453145 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453146

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/

echo "* Start cufflinks [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453146 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453147

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/

echo "* Start cufflinks [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453147 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453148

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/

echo "* Start cufflinks [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453148 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453149

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/

echo "* Start cufflinks [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453149 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


