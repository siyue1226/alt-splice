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
# lane SRR070411

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/

echo "* Start tophat [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/trimmed/SRR070411_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070428

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/

echo "* Start tophat [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/trimmed/SRR070428_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070411

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/

echo "* Start cufflinks [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070411 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_fat_body] [SRR070411] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070428

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/

echo "* Start cufflinks [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_fat_body SRR070428 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_fat_body] [SRR070428] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


