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
# lane SRR453099

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/

echo "* Start tophat [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/trimmed/SRR453099_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453100

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/

echo "* Start tophat [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/trimmed/SRR453100_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453101

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/

echo "* Start tophat [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/trimmed/SRR453101_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453102

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/

echo "* Start tophat [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/trimmed/SRR453102_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453103

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/

echo "* Start tophat [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/trimmed/SRR453103_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453104

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/

echo "* Start tophat [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/trimmed/SRR453104_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453105

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/

echo "* Start tophat [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/trimmed/SRR453105_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453106

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/

echo "* Start tophat [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/trimmed/SRR453106_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453107

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/

echo "* Start tophat [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/trimmed/SRR453107_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453108

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/

echo "* Start tophat [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/trimmed/SRR453108_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453099

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/

echo "* Start cufflinks [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453099 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453099] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453100

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/

echo "* Start cufflinks [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453100 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453100] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453101

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/

echo "* Start cufflinks [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453101 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453101] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453102

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/

echo "* Start cufflinks [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453102 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453102] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453103

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/

echo "* Start cufflinks [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453103 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453103] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453104

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/

echo "* Start cufflinks [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453104 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453104] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453105

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/

echo "* Start cufflinks [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453105 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453105] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453106

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/

echo "* Start cufflinks [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453106 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453106] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453107

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/

echo "* Start cufflinks [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453107 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453107] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453108

cd /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/

echo "* Start cufflinks [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` SmIntestine SRR453108 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [SmIntestine] [SRR453108] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


