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
# lane SRR070397

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/

echo "* Start tophat [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/trimmed/SRR070397_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR182358

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/

echo "* Start tophat [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/trimmed/SRR182358_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR350959

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/

echo "* Start tophat [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/trimmed/SRR350959_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070397

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/

echo "* Start cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR070397 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR070397] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR182358

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/

echo "* Start cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR182358 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR182358] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR350959

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/

echo "* Start cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_accessory_glands SRR350959 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_accessory_glands] [SRR350959] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


