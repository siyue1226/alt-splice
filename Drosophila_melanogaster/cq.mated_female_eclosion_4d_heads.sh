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
# lane SRR070414

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/

echo "* Start cuffquant [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070415

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/

echo "* Start cuffquant [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


