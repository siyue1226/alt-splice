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
# lane SRR070421

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/

echo "* Start tophat [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/trimmed/SRR070421_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070424

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/

echo "* Start tophat [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/trimmed/SRR070424_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070421

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/

echo "* Start cufflinks [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070421 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_20d_heads] [SRR070421] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070424

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/

echo "* Start cufflinks [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_20d_heads SRR070424 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_20d_heads] [SRR070424] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


