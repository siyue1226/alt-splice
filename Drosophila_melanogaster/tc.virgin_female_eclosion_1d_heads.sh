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
# lane SRR070436

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/

echo "* Start tophat [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/trimmed/SRR070436_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070437

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/

echo "* Start tophat [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/trimmed/SRR070437_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100281

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/

echo "* Start tophat [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/trimmed/SRR100281_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070436

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/

echo "* Start cufflinks [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070436 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_1d_heads] [SRR070436] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070437

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/

echo "* Start cufflinks [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR070437 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_1d_heads] [SRR070437] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100281

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/

echo "* Start cufflinks [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_1d_heads SRR100281 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_1d_heads] [SRR100281] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


