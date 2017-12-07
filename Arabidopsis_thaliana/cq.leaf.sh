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

#----------------------------#
# cuffquant
#----------------------------#
# lane SRR6179904

cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179904/

echo "* Start cuffquant [leaf] [SRR6179904] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.mask.gtf \
    -b /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa \
    /home/wangq/data/rna-seq/Ath_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/Ath_trans/process/SRR6179904/cq_out \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179904/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` leaf SRR6179904 [cuffquant] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End cuffquant [leaf] [SRR6179904] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffquant.log
