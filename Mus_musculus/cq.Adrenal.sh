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
# lane SRR453116

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/

echo "* Start cuffquant [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453117

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/

echo "* Start cuffquant [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453118

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/

echo "* Start cuffquant [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453119

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/

echo "* Start cuffquant [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453120

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/

echo "* Start cuffquant [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453121

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/

echo "* Start cuffquant [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


