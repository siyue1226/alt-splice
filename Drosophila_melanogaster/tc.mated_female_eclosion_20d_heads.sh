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
# lane SRR070420

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/

echo "* Start tophat [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/trimmed/SRR070420_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100274

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/

echo "* Start tophat [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/trimmed/SRR100274_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR116383

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/

echo "* Start tophat [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/trimmed/SRR116383_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111882

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/

echo "* Start tophat [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/trimmed/SRR111882_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070420

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/

echo "* Start cufflinks [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR070420 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_20d_heads] [SRR070420] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100274

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/

echo "* Start cufflinks [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR100274 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_20d_heads] [SRR100274] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR116383

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/

echo "* Start cufflinks [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR116383 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_20d_heads] [SRR116383] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111882

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/

echo "* Start cufflinks [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_female_eclosion_20d_heads SRR111882 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_female_eclosion_20d_heads] [SRR111882] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


