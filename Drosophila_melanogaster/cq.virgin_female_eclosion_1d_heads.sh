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
# lane SRR070436

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/

echo "* Start cuffquant [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070437

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/

echo "* Start cuffquant [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100281

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/

echo "* Start cuffquant [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


