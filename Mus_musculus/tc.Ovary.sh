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
# lane SRR453077

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/

echo "* Start tophat [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/trimmed/SRR453077_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453077 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453078

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/

echo "* Start tophat [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/trimmed/SRR453078_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453078 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453079

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/

echo "* Start tophat [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/trimmed/SRR453079_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453079 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453080

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/

echo "* Start tophat [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/trimmed/SRR453080_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453080 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453081

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/

echo "* Start tophat [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/trimmed/SRR453081_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453081 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453082

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/

echo "* Start tophat [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/trimmed/SRR453082_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453082 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453083

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/

echo "* Start tophat [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/trimmed/SRR453083_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453083 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453084

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/

echo "* Start tophat [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/trimmed/SRR453084_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453084 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453085

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/

echo "* Start tophat [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/trimmed/SRR453085_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453085 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453086

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/

echo "* Start tophat [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/trimmed/SRR453086_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Ovary SRR453086 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453077

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/

echo "* Start cufflinks [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453077 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453077] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453078

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/

echo "* Start cufflinks [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453078 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453078] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453079

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/

echo "* Start cufflinks [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453079 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453079] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453080

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/

echo "* Start cufflinks [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453080 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453080] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453081

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/

echo "* Start cufflinks [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453081 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453081] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453082

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/

echo "* Start cufflinks [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453082 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453082] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453083

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/

echo "* Start cufflinks [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453083 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453083] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453084

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/

echo "* Start cufflinks [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453084 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453084] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453085

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/

echo "* Start cufflinks [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453085 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453085] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453086

cd /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/

echo "* Start cufflinks [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Ovary SRR453086 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Ovary] [SRR453086] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


