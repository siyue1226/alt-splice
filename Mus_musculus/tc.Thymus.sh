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
# tophat
#----------------------------#
# lane SRR453134

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/

echo "* Start tophat [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/trimmed/SRR453134_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453134 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453135

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/

echo "* Start tophat [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/trimmed/SRR453135_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453135 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453136

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/

echo "* Start tophat [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/trimmed/SRR453136_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453136 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453137

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/

echo "* Start tophat [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/trimmed/SRR453137_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453137 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453138

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/

echo "* Start tophat [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/trimmed/SRR453138_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453138 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453139

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/

echo "* Start tophat [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/trimmed/SRR453139_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Thymus SRR453139 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453134

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/

echo "* Start cufflinks [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453134 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453134] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453135

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/

echo "* Start cufflinks [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453135 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453135] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453136

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/

echo "* Start cufflinks [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453136 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453136] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453137

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/

echo "* Start cufflinks [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453137 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453137] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453138

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/

echo "* Start cufflinks [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453138 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453138] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453139

cd /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/

echo "* Start cufflinks [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Thymus SRR453139 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Thymus] [SRR453139] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


