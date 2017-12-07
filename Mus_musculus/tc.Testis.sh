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
# lane SRR453140

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/

echo "* Start tophat [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/trimmed/SRR453140_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453140 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453141

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/

echo "* Start tophat [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/trimmed/SRR453141_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453141 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453142

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/

echo "* Start tophat [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/trimmed/SRR453142_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453142 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# lane SRR453143

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/

echo "* Start tophat [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/th_out \
    /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65 \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_1.sickle.fq.gz \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/trimmed/SRR453143_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` Testis SRR453143 [tophat] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End tophat [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane SRR453140

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/

echo "* Start cufflinks [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453140 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Testis] [SRR453140] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453141

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/

echo "* Start cufflinks [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453141 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Testis] [SRR453141] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453142

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/

echo "* Start cufflinks [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453142 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Testis] [SRR453142] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# lane SRR453143

cd /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/

echo "* Start cufflinks [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check -M /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.mask.gtf \
    -o /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/cl_out \
    /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` Testis SRR453143 [cufflinks] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cufflinks [Testis] [SRR453143] `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cufflinks.log


