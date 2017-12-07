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
# lane SRR070396

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/

echo "* Start tophat [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/trimmed/SRR070396_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070417

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/

echo "* Start tophat [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/trimmed/SRR070417_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070396

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/

echo "* Start cufflinks [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070396 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_4d_ovaries] [SRR070396] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070417

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/

echo "* Start cufflinks [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_ovaries SRR070417 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [virgin_female_eclosion_4d_ovaries] [SRR070417] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


