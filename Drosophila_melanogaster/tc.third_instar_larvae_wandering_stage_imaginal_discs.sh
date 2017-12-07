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
# tophat
#----------------------------#
# lane SRR070392

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/trimmed/SRR070392_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070393

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/trimmed/SRR070393_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111884

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/trimmed/SRR111884_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR111885

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/trimmed/SRR111885_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR350962

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/trimmed/SRR350962_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR350963

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/

echo "* Start tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/trimmed/SRR350963_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070392

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070392 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070392] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070393

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR070393 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR070393] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111884

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111884 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111884] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR111885

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR111885 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR111885] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR350962

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350962 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350962] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR350963

cd /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/

echo "* Start cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` third_instar_larvae_wandering_stage_imaginal_discs SRR350963 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [third_instar_larvae_wandering_stage_imaginal_discs] [SRR350963] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


