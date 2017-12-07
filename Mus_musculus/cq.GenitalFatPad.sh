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
# lane SRR453126

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/

echo "* Start cuffquant [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453126 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [GenitalFatPad] [SRR453126] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453127

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/

echo "* Start cuffquant [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453127 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [GenitalFatPad] [SRR453127] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453128

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/

echo "* Start cuffquant [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453128 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [GenitalFatPad] [SRR453128] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453129

cd /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/

echo "* Start cuffquant [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` GenitalFatPad SRR453129 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [GenitalFatPad] [SRR453129] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


