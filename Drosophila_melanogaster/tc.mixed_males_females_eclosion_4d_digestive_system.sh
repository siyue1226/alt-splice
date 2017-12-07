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
# lane SRR070401

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/

echo "* Start tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/trimmed/SRR070401_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111878

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/

echo "* Start tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/trimmed/SRR111878_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111879

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/

echo "* Start tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/trimmed/SRR111879_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070401

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/

echo "* Start cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR070401 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR070401] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111878

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/

echo "* Start cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111878 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR111878] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111879

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/

echo "* Start cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_4d_digestive_system SRR111879 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_4d_digestive_system] [SRR111879] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


