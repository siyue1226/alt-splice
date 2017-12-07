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
# lane SRR070414

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/

echo "* Start tophat [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/trimmed/SRR070414_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070415

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/

echo "* Start tophat [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/trimmed/SRR070415_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070414

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/

echo "* Start cufflinks [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070414 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_4d_heads] [SRR070414] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070415

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/

echo "* Start cufflinks [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_4d_heads SRR070415 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_4d_heads] [SRR070415] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


