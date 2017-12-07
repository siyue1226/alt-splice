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
# lane SRR453099

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/

echo "* Start cuffquant [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453100

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/

echo "* Start cuffquant [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453101

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/

echo "* Start cuffquant [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453102

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/

echo "* Start cuffquant [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453103

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/

echo "* Start cuffquant [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453104

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/

echo "* Start cuffquant [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453105

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/

echo "* Start cuffquant [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453106

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/

echo "* Start cuffquant [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453107

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/

echo "* Start cuffquant [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453108

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/

echo "* Start cuffquant [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


