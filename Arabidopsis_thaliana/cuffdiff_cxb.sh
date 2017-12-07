#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/Ath_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29

cd /home/wangq/data/rna-seq/Ath_trans/process

#----------------------------#
# cuffdiff_cxb
#----------------------------#
echo "* Start cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffmergediff.log

/home/wangq/bin/cuffdiff -p 14 \
    --no-update-check -u -M /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.mask.gtf \
    -b /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa \
    /home/wangq/data/rna-seq/Ath_trans/process/merged_asm/merged.gtf \
    -L flower,root,leaf,steam,seedling \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179906/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179904/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/cq_out/abundances.cxb \
    -o /home/wangq/data/rna-seq/Ath_trans/process/diff_out_cxb

[ $? -ne 0 ] && echo `date` ] [cuffdiff_cxb] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffmergediff.log

