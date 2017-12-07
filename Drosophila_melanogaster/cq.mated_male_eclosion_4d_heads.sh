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
# lane SRR070400

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/

echo "* Start cuffquant [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070416

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/

echo "* Start cuffquant [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


