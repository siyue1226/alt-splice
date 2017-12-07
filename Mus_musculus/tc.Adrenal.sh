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
# lane SRR453116

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/

echo "* Start tophat [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/trimmed/SRR453116_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453117

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/

echo "* Start tophat [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/trimmed/SRR453117_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453118

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/

echo "* Start tophat [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/trimmed/SRR453118_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453119

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/

echo "* Start tophat [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/trimmed/SRR453119_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453120

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/

echo "* Start tophat [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/trimmed/SRR453120_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453121

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/

echo "* Start tophat [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/trimmed/SRR453121_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453116

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/

echo "* Start cufflinks [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453116 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453116] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453117

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/

echo "* Start cufflinks [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453117 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453117] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453118

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/

echo "* Start cufflinks [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453118 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453118] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453119

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/

echo "* Start cufflinks [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453119 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453119] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453120

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/

echo "* Start cufflinks [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453120 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453120] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453121

cd /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/

echo "* Start cufflinks [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Adrenal SRR453121 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Adrenal] [SRR453121] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


