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
# lane SRR070399

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/

echo "* Start cuffquant [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR070395

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/

echo "* Start cuffquant [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


