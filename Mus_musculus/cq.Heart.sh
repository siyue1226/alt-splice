#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/mouse_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65

#----------------------------#
# cuffquant
#----------------------------#
# lane SRR453172

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/

echo "* Start cuffquant [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453172 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Heart] [SRR453172] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453173

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/

echo "* Start cuffquant [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453173 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Heart] [SRR453173] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453174

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/

echo "* Start cuffquant [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453174 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Heart] [SRR453174] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453175

cd /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/

echo "* Start cuffquant [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Heart SRR453175 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Heart] [SRR453175] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


