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
# lane SRR070404

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/

echo "* Start tophat [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/trimmed/SRR070404_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070391

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/

echo "* Start tophat [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/trimmed/SRR070391_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070404

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/

echo "* Start cufflinks [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070404 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_20d_carcass] [SRR070404] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070391

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/

echo "* Start cufflinks [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_carcass SRR070391 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_20d_carcass] [SRR070391] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


