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
# tophat
#----------------------------#
# lane SRR070408

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/

echo "* Start tophat [third_instar_larvae_wandering_stage_digestive_system] [SRR070408] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/trimmed/SRR070408_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/trimmed/SRR070408_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_digestive_system SRR070408 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_digestive_system] [SRR070408] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100268

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/

echo "* Start tophat [third_instar_larvae_wandering_stage_digestive_system] [SRR100268] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/trimmed/SRR100268_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/trimmed/SRR100268_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_digestive_system SRR100268 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_digestive_system] [SRR100268] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070408

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_digestive_system] [SRR070408] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_digestive_system SRR070408 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_digestive_system] [SRR070408] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100268

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_digestive_system] [SRR100268] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_digestive_system SRR100268 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_digestive_system] [SRR100268] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


