#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/rice_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa /home/wangq/data/rna-seq/rice_trans/ref/rice.29

cd /home/wangq/data/rna-seq/rice_trans/process

#----------------------------#
# cuffdiff_cxb
#----------------------------#
echo "* Start cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cuffmergediff.log

/home/wangq/bin/cuffdiff -p 14 \
    --no-update-check -u -M /home/wangq/data/rna-seq/rice_trans/ref/rice.29.mask.gtf \
    -b /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa \
    /home/wangq/data/rna-seq/rice_trans/process/merged_asm/merged.gtf \
    -L AA,AB,AC,AD,AE,AF,AG,AH,AK,BH,CA \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352184/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352187/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352189/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352190/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352192/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352194/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352204/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352206/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352207/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352209/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/rice_trans/process/SRR352211/cq_out/abundances.cxb \
    

    -o /home/wangq/data/rna-seq/rice_trans/process/diff_out_cxb

[ $? -ne 0 ] && echo `date` ] [cuffdiff_cxb] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cuffmergediff.log

