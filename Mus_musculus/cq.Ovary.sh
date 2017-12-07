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
# lane SRR453077

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/

echo "* Start cuffquant [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453077 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453078

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/

echo "* Start cuffquant [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453078 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453079

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/

echo "* Start cuffquant [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453079 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453080

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/

echo "* Start cuffquant [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453080 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453081

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/

echo "* Start cuffquant [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453081 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453082

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/

echo "* Start cuffquant [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453082 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453083

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/

echo "* Start cuffquant [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453083 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453084

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/

echo "* Start cuffquant [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453084 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453085

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/

echo "* Start cuffquant [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453085 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# lane SRR453086

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/

echo "* Start cuffquant [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log

# cuffquant
/home/wangq/bin/cuffquant -p 8 \
    --no-update-check -u -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -b /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    /home/wangq/data/rna-seq/mouse_trans/process/merged_asm/merged.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/cq_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453086 [cuffquant] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffquant [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffquant.log


