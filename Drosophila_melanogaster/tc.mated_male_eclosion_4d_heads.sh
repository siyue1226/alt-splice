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
# lane SRR070400

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/

echo "* Start tophat [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/trimmed/SRR070400_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070416

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/

echo "* Start tophat [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/trimmed/SRR070416_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070400

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/

echo "* Start cufflinks [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070400 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_heads] [SRR070400] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070416

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/

echo "* Start cufflinks [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_heads SRR070416 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_heads] [SRR070416] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


