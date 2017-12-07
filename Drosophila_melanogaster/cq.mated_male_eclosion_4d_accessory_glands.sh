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
# lane SRR070397

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/

echo "* Start cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR182358

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/

echo "* Start cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR350959

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/

echo "* Start cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


