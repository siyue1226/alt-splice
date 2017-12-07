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
# lane SRR453087

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/

echo "* Start tophat [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/trimmed/SRR453087_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453088

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/

echo "* Start tophat [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/trimmed/SRR453088_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453089

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/

echo "* Start tophat [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/trimmed/SRR453089_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453090

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/

echo "* Start tophat [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/trimmed/SRR453090_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453091

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/

echo "* Start tophat [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/trimmed/SRR453091_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453092

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/

echo "* Start tophat [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/trimmed/SRR453092_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453087

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/

echo "* Start cufflinks [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453087 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453087] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453088

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/

echo "* Start cufflinks [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453088 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453088] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453089

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/

echo "* Start cufflinks [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453089 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453089] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453090

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/

echo "* Start cufflinks [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453090 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453090] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453091

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/

echo "* Start cufflinks [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453091 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453091] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453092

cd /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/

echo "* Start cufflinks [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` MammaryGland SRR453092 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [MammaryGland] [SRR453092] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


