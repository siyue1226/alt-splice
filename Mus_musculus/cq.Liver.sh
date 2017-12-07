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
# lane SRR453150

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/

echo "* Start cuffquant [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453150 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453150] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453151

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/

echo "* Start cuffquant [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453151 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453151] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453152

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/

echo "* Start cuffquant [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453152 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453152] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453153

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/

echo "* Start cuffquant [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453153 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453153] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453154

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/

echo "* Start cuffquant [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453154 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453154] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453155

cd /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/

echo "* Start cuffquant [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Liver SRR453155 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Liver] [SRR453155] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


