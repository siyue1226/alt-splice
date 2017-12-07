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
# lane SRR070422

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/

echo "* Start tophat [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/trimmed/SRR070422_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR070423

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/

echo "* Start tophat [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/trimmed/SRR070423_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR100276

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/

echo "* Start tophat [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/trimmed/SRR100276_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR350960

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/

echo "* Start tophat [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/trimmed/SRR350960_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# lane SRR350961

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/

echo "* Start tophat [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/th_out \
    /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65 \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/trimmed/SRR350961_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [tophat] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End tophat [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR070422

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/

echo "* Start cufflinks [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070422 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_testes] [SRR070422] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR070423

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/

echo "* Start cufflinks [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR070423 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_testes] [SRR070423] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR100276

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/

echo "* Start cufflinks [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR100276 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_testes] [SRR100276] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR350960

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/

echo "* Start cufflinks [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350960 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_testes] [SRR350960] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# lane SRR350961

cd /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/

echo "* Start cufflinks [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -o /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/cl_out \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` mated_male_eclosion_4d_testes SRR350961 [cufflinks] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cufflinks [mated_male_eclosion_4d_testes] [SRR350961] `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cufflinks.log


