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
# lane SRR070427

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/

echo "* Start tophat [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/trimmed/SRR070427_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/trimmed/SRR070427_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR070427 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100270

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/

echo "* Start tophat [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/trimmed/SRR100270_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/trimmed/SRR100270_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR100270 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070427

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/

echo "* Start cufflinks [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR070427 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_salivary_glands] [SRR070427] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100270

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/

echo "* Start cufflinks [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_salivary_glands SRR100270 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_salivary_glands] [SRR100270] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


