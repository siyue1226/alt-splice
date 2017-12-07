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
# lane SRR070388

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/

echo "* Start tophat [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/trimmed/SRR070388_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070419

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/

echo "* Start tophat [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/trimmed/SRR070419_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100275

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/

echo "* Start tophat [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/trimmed/SRR100275_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070388

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/

echo "* Start cufflinks [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070388 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_20d_heads] [SRR070388] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070419

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/

echo "* Start cufflinks [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR070419 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_20d_heads] [SRR070419] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100275

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/

echo "* Start cufflinks [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_20d_heads SRR100275 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_20d_heads] [SRR100275] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


