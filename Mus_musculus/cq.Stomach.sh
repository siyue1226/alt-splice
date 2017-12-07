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
# lane SRR453093

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/

echo "* Start cuffquant [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453093 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453093] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453094

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/

echo "* Start cuffquant [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453094 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453094] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453095

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/

echo "* Start cuffquant [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453095 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453095] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453096

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/

echo "* Start cuffquant [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453096 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453096] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453097

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/

echo "* Start cuffquant [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453097 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453097] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453098

cd /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/

echo "* Start cuffquant [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Stomach SRR453098 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Stomach] [SRR453098] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


