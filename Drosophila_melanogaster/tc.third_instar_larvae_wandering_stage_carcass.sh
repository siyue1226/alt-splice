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
# lane SRR070426

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/

echo "* Start tophat [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/trimmed/SRR070426_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100269

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/

echo "* Start tophat [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/trimmed/SRR100269_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070426

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR070426 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_carcass] [SRR070426] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100269

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_carcass SRR100269 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_carcass] [SRR100269] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


