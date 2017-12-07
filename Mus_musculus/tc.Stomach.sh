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
# lane SRR453093

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/

echo "* Start tophat [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/trimmed/SRR453093_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453093 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453094

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/

echo "* Start tophat [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/trimmed/SRR453094_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453094 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453095

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/

echo "* Start tophat [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/trimmed/SRR453095_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453095 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453096

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/

echo "* Start tophat [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/trimmed/SRR453096_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453096 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453097

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/

echo "* Start tophat [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/trimmed/SRR453097_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453097 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453098

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/

echo "* Start tophat [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/trimmed/SRR453098_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Stomach SRR453098 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453093

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/

echo "* Start cufflinks [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453093 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453094

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/

echo "* Start cufflinks [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453094 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453095

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/

echo "* Start cufflinks [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453095 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453096

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/

echo "* Start cufflinks [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453096 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453097

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/

echo "* Start cufflinks [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453097 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453098

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/

echo "* Start cufflinks [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453098 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


