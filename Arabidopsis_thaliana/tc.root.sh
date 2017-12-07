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
# tophat
#----------------------------#
# lane SRR6179907

cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/

echo "* Start tophat [root] [SRR6179907] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.gtf \
    -o /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/th_out \
    /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29 \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/trimmed/SRR6179907_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/trimmed/SRR6179907_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` root SRR6179907 [tophat] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End tophat [root] [SRR6179907] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/tophat.log




#----------------------------#
# cufflinks
#----------------------------#
# lane SRR6179907

cd /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/

echo "* Start cufflinks [root] [SRR6179907] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.mask.gtf \
    -o /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/cl_out \
    /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` root SRR6179907 [cufflinks] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End cufflinks [root] [SRR6179907] `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cufflinks.log

