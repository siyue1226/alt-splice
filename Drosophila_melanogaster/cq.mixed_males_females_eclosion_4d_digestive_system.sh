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
# lane SRR070401

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/

echo "* Start cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR111878

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/

echo "* Start cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR111879

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/

echo "* Start cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


