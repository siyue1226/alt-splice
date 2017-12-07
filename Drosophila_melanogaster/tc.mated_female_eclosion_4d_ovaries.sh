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
# lane SRR070431

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/

echo "* Start tophat [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/trimmed/SRR070431_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100277

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/

echo "* Start tophat [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/trimmed/SRR100277_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100283

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/

echo "* Start tophat [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/trimmed/SRR100283_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070431

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/

echo "* Start cufflinks [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR070431 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_4d_ovaries] [SRR070431] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100277

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/

echo "* Start cufflinks [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100277 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_4d_ovaries] [SRR100277] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100283

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/

echo "* Start cufflinks [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_ovaries SRR100283 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_4d_ovaries] [SRR100283] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


