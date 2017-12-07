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
# lane SRR070434

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/

echo "* Start cuffquant [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070434 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_1d_heads] [SRR070434] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070435

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/

echo "* Start cuffquant [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR070435 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_1d_heads] [SRR070435] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100279

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/

echo "* Start cuffquant [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_1d_heads SRR100279 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_1d_heads] [SRR100279] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


