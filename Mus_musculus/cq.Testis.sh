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
# lane SRR453140

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/

echo "* Start cuffquant [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453140 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453141

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/

echo "* Start cuffquant [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453141 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453142

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/

echo "* Start cuffquant [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453142 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453143

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/

echo "* Start cuffquant [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453143 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


