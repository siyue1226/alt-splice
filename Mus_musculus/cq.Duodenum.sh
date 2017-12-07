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
# cuffquant
#----------------------------#
# lane SRR453109

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/

echo "* Start cuffquant [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453109 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453109] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453110

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/

echo "* Start cuffquant [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453110 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453110] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453111

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/

echo "* Start cuffquant [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453111 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453111] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453112

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/

echo "* Start cuffquant [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453112 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453112] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453113

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/

echo "* Start cuffquant [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453113 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453113] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453114

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/

echo "* Start cuffquant [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453114 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453114] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453115

cd /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/

echo "* Start cuffquant [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Duodenum SRR453115 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Duodenum] [SRR453115] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


