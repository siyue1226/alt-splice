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
# lane SRR070403

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/

echo "* Start tophat [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/trimmed/SRR070403_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111883

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/

echo "* Start tophat [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/trimmed/SRR111883_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070403

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/

echo "* Start cufflinks [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR070403 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_20d_digestive_system] [SRR070403] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111883

cd /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/

echo "* Start cufflinks [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mixed_males_females_eclosion_20d_digestive_system SRR111883 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mixed_males_females_eclosion_20d_digestive_system] [SRR111883] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


