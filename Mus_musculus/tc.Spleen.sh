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
# lane SRR453160

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/

echo "* Start tophat [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/trimmed/SRR453160_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453160 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453161

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/

echo "* Start tophat [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/trimmed/SRR453161_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453161 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453162

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/

echo "* Start tophat [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/trimmed/SRR453162_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453162 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453163

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/

echo "* Start tophat [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/trimmed/SRR453163_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453163 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453164

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/

echo "* Start tophat [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/trimmed/SRR453164_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453164 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453165

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/

echo "* Start tophat [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/trimmed/SRR453165_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Spleen SRR453165 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453160

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/

echo "* Start cufflinks [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453160 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453160] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453161

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/

echo "* Start cufflinks [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453161 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453161] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453162

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/

echo "* Start cufflinks [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453162 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453162] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453163

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/

echo "* Start cufflinks [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453163 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453163] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453164

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/

echo "* Start cufflinks [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453164 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453164] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453165

cd /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/

echo "* Start cufflinks [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Spleen SRR453165 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Spleen] [SRR453165] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


