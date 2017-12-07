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
# lane SRR453156

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/

echo "* Start cuffquant [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453156 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Lung] [SRR453156] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453157

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/

echo "* Start cuffquant [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453157 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Lung] [SRR453157] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453158

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/

echo "* Start cuffquant [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453158 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Lung] [SRR453158] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453159

cd /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/

echo "* Start cuffquant [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Lung SRR453159 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Lung] [SRR453159] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


