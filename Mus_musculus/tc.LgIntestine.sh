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
# lane SRR453122

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/

echo "* Start tophat [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/trimmed/SRR453122_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453123

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/

echo "* Start tophat [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/trimmed/SRR453123_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453124

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/

echo "* Start tophat [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/trimmed/SRR453124_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453125

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/

echo "* Start tophat [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/trimmed/SRR453125_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453122

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/

echo "* Start cufflinks [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453123

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/

echo "* Start cufflinks [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453124

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/

echo "* Start cufflinks [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453125

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/

echo "* Start cufflinks [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


