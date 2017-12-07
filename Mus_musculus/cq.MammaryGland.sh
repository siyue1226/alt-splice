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
# lane SRR453087

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/

echo "* Start cuffquant [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453088

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/

echo "* Start cuffquant [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453089

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/

echo "* Start cuffquant [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453090

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/

echo "* Start cuffquant [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453091

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/

echo "* Start cuffquant [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453092

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/

echo "* Start cuffquant [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


