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
# lane SRR453150

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/

echo "* Start tophat [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/trimmed/SRR453150_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453150 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453151

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/

echo "* Start tophat [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/trimmed/SRR453151_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453151 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453152

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/

echo "* Start tophat [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/trimmed/SRR453152_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453152 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453153

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/

echo "* Start tophat [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/trimmed/SRR453153_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453153 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453154

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/

echo "* Start tophat [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/trimmed/SRR453154_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453154 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453155

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/

echo "* Start tophat [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/trimmed/SRR453155_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Liver SRR453155 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453150

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/

echo "* Start cufflinks [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453150 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453151

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/

echo "* Start cufflinks [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453151 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453152

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/

echo "* Start cufflinks [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453152 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453153

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/

echo "* Start cufflinks [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453153 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453154

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/

echo "* Start cufflinks [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453154 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453155

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/

echo "* Start cufflinks [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453155 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


