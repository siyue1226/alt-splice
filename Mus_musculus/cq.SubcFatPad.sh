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
# lane SRR453130

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/

echo "* Start cuffquant [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453131

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/

echo "* Start cuffquant [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453132

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/

echo "* Start cuffquant [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453133

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/

echo "* Start cuffquant [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


