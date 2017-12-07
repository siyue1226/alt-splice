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
# lane SRR070407

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/

echo "* Start tophat [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/trimmed/SRR070407_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070425

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/

echo "* Start tophat [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/trimmed/SRR070425_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070407

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070407 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_salivary_glands] [SRR070407] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070425

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_salivary_glands SRR070425 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_salivary_glands] [SRR070425] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


