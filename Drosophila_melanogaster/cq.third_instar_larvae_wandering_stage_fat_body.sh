#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/dmel_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65

#----------------------------#
# cuffquant
#----------------------------#
# lane SRR070405

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070405/

echo "* Start cuffquant [third_instar_larvae_wandering_stage_fat_body] [SRR070405] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070405/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070405/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_fat_body SRR070405 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [third_instar_larvae_wandering_stage_fat_body] [SRR070405] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070406

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070406/

echo "* Start cuffquant [third_instar_larvae_wandering_stage_fat_body] [SRR070406] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070406/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070406/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_fat_body SRR070406 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [third_instar_larvae_wandering_stage_fat_body] [SRR070406] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


