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
# lane SRR070412

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/

echo "* Start tophat [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/trimmed/SRR070412_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100271

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/

echo "* Start tophat [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/trimmed/SRR100271_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070412

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/

echo "* Start cufflinks [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR070412 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_2d_CNS] [SRR070412] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100271

cd /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/

echo "* Start cufflinks [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` WPP_2d_CNS SRR100271 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [WPP_2d_CNS] [SRR100271] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


