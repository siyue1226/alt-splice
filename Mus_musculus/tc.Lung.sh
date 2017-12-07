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
# lane SRR453156

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/

echo "* Start tophat [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/trimmed/SRR453156_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453156 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453157

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/

echo "* Start tophat [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/trimmed/SRR453157_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453157 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453158

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/

echo "* Start tophat [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/trimmed/SRR453158_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453158 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453159

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/

echo "* Start tophat [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/trimmed/SRR453159_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Lung SRR453159 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453156

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/

echo "* Start cufflinks [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453156 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453157

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/

echo "* Start cufflinks [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453157 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453158

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/

echo "* Start cufflinks [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453158 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453159

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/

echo "* Start cufflinks [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453159 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


