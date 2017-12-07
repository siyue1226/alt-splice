#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/rice_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa /home/wangq/data/rna-seq/rice_trans/ref/rice.29

#----------------------------#
# tophat
#----------------------------#
# lane SRR352194

cd /home/wangq/data/rna-seq/rice_trans/process/SRR352194/

echo "* Start tophat [OSN_AF] [SRR352194] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/rice_trans/ref/rice.29.gtf \
    -o /home/wangq/data/rna-seq/rice_trans/process/SRR352194/th_out \
    /home/wangq/data/rna-seq/rice_trans/ref/rice.29 \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352194/trimmed/SRR352194_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352194/trimmed/SRR352194_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` OSN_AF SRR352194 [tophat] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End tophat [OSN_AF] [SRR352194] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/tophat.log




#----------------------------#
# cufflinks
#----------------------------#
# lane SRR352194

cd /home/wangq/data/rna-seq/rice_trans/process/SRR352194/

echo "* Start cufflinks [OSN_AF] [SRR352194] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/rice_trans/ref/rice.29.mask.gtf \
    -o /home/wangq/data/rna-seq/rice_trans/process/SRR352194/cl_out \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352194/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` OSN_AF SRR352194 [cufflinks] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End cufflinks [OSN_AF] [SRR352194] `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cufflinks.log

