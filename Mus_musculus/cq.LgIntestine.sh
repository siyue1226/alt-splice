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
# lane SRR453122

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/

echo "* Start cuffquant [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453122 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [LgIntestine] [SRR453122] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453123

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/

echo "* Start cuffquant [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453123 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [LgIntestine] [SRR453123] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453124

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/

echo "* Start cuffquant [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453124 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [LgIntestine] [SRR453124] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453125

cd /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/

echo "* Start cuffquant [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` LgIntestine SRR453125 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [LgIntestine] [SRR453125] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


