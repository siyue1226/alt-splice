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
# lane SRR453166

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/

echo "* Start cuffquant [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453166 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453166] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453167

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/

echo "* Start cuffquant [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453167 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453167] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453168

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/

echo "* Start cuffquant [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453168 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453168] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453169

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/

echo "* Start cuffquant [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453169 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453169] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453170

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/

echo "* Start cuffquant [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453170 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453170] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453171

cd /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/

echo "* Start cuffquant [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Colon SRR453171 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Colon] [SRR453171] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


