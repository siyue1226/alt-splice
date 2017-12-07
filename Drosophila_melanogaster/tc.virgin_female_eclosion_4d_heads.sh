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
# lane SRR070430

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/

echo "* Start tophat [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/trimmed/SRR070430_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100278

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/

echo "* Start tophat [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/trimmed/SRR100278_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100282

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/

echo "* Start tophat [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/trimmed/SRR100282_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070430

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/

echo "* Start cufflinks [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100278

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/

echo "* Start cufflinks [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100282

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/

echo "* Start cufflinks [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


