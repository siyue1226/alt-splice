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
# lane SRR453172

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/

echo "* Start tophat [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/trimmed/SRR453172_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453172 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453173

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/

echo "* Start tophat [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/trimmed/SRR453173_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453173 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453174

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/

echo "* Start tophat [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/trimmed/SRR453174_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453174 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453175

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/

echo "* Start tophat [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/trimmed/SRR453175_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Heart SRR453175 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453172

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/

echo "* Start cufflinks [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453172 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453173

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/

echo "* Start cufflinks [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453173 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453174

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/

echo "* Start cufflinks [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453174 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453175

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/

echo "* Start cufflinks [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453175 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


