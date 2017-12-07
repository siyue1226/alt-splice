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
# lane SRR070427

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/

echo "* Start cuffquant [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR070427 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100270

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/

echo "* Start cuffquant [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR100270 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


