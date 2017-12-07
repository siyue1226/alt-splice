#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/dmel_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65

#----------------------------#
# cuffquant
#----------------------------#
# lane SRR070430

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/

echo "* Start cuffquant [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR070430 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_4d_heads] [SRR070430] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100278

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/

echo "* Start cuffquant [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100278 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_4d_heads] [SRR100278] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# lane SRR100282

cd /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/

echo "* Start cuffquant [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/cq_out \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` virgin_female_eclosion_4d_heads SRR100282 [cuffquant] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffquant [virgin_female_eclosion_4d_heads] [SRR100282] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffquant.log


