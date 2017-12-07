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
# lane SRR453134

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/

echo "* Start cuffquant [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453134 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453135

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/

echo "* Start cuffquant [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453135 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453136

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/

echo "* Start cuffquant [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453136 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453137

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/

echo "* Start cuffquant [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453137 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453138

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/

echo "* Start cuffquant [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453138 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453139

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/

echo "* Start cuffquant [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453139 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


