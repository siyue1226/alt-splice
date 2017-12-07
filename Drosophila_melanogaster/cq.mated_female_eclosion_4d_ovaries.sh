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
# lane SRR070431

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/

echo "* Start cuffquant [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100277

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/

echo "* Start cuffquant [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100283

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/

echo "* Start cuffquant [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


