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
# lane SRR453109

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/

echo "* Start tophat [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/trimmed/SRR453109_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453110

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/

echo "* Start tophat [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/trimmed/SRR453110_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453111

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/

echo "* Start tophat [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/trimmed/SRR453111_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453112

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/

echo "* Start tophat [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/trimmed/SRR453112_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453113

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/

echo "* Start tophat [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/trimmed/SRR453113_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453114

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/

echo "* Start tophat [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/trimmed/SRR453114_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453115

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/

echo "* Start tophat [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/trimmed/SRR453115_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453109

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/

echo "* Start cufflinks [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453110

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/

echo "* Start cufflinks [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453111

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/

echo "* Start cufflinks [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453112

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/

echo "* Start cufflinks [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453113

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/

echo "* Start cufflinks [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453114

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/

echo "* Start cufflinks [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453115

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/

echo "* Start cufflinks [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


