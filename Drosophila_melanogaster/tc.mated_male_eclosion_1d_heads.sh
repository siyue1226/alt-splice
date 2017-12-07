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
# lane SRR070432

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/

echo "* Start tophat [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/trimmed/SRR070432_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070433

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/

echo "* Start tophat [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/trimmed/SRR070433_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100280

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/

echo "* Start tophat [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/trimmed/SRR100280_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070432

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/

echo "* Start cufflinks [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070432 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_1d_heads] [SRR070432] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070433

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/

echo "* Start cufflinks [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR070433 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_1d_heads] [SRR070433] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100280

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/

echo "* Start cufflinks [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_1d_heads SRR100280 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_1d_heads] [SRR100280] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


