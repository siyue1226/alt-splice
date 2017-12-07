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
# lane SRR453144

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/

echo "* Start cuffquant [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453144 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453144] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453145

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/

echo "* Start cuffquant [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453145 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453145] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453146

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/

echo "* Start cuffquant [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453146 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453146] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453147

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/

echo "* Start cuffquant [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453147 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453147] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453148

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/

echo "* Start cuffquant [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453148 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453148] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453149

cd /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/

echo "* Start cuffquant [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Kidney SRR453149 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Kidney] [SRR453149] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


