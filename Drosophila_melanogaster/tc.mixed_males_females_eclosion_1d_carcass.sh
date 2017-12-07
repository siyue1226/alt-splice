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
# lane SRR070399

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/

echo "* Start tophat [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/trimmed/SRR070399_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070395

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/

echo "* Start tophat [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/trimmed/SRR070395_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070399

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/

echo "* Start cufflinks [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070399 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_1d_carcass] [SRR070399] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070395

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/

echo "* Start cufflinks [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_1d_carcass SRR070395 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_1d_carcass] [SRR070395] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


