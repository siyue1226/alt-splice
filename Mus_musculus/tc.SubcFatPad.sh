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
# lane SRR453130

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/

echo "* Start tophat [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/trimmed/SRR453130_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453131

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/

echo "* Start tophat [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/trimmed/SRR453131_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453132

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/

echo "* Start tophat [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/trimmed/SRR453132_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453133

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/

echo "* Start tophat [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/trimmed/SRR453133_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453130

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/

echo "* Start cufflinks [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453130 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SubcFatPad] [SRR453130] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453131

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/

echo "* Start cufflinks [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453131 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SubcFatPad] [SRR453131] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453132

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/

echo "* Start cufflinks [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453132 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SubcFatPad] [SRR453132] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453133

cd /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/

echo "* Start cufflinks [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SubcFatPad SRR453133 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SubcFatPad] [SRR453133] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


