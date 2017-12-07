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
# lane SRR070412

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/

echo "* Start cuffquant [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100271

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/

echo "* Start cuffquant [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


