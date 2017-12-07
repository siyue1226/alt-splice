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
# lane SRR453160

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/

echo "* Start cuffquant [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453160 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453161

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/

echo "* Start cuffquant [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453161 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453162

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/

echo "* Start cuffquant [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453162 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453163

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/

echo "* Start cuffquant [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453163 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453164

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/

echo "* Start cuffquant [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453164 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453165

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/

echo "* Start cuffquant [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453165 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


