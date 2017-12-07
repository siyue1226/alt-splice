#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/bodymap2

### bowtie index
# bowtie2-build /home/wangq/data/bodymap2/ref/human.37.fa /home/wangq/data/bodymap2/ref/human.37

#----------------------------#
# tophat
#----------------------------#
# lane ERR030877

cd /home/wangq/data/bodymap2/process/prostate/ERR030877/

echo "* Start tophat [prostate] [ERR030877] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/prostate/ERR030877/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/prostate/ERR030877/trimmed/ERR030877_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` prostate ERR030877 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [prostate] [ERR030877] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030898

cd /home/wangq/data/bodymap2/process/prostate/ERR030898/

echo "* Start tophat [prostate] [ERR030898] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/prostate/ERR030898/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/prostate/ERR030898/trimmed/ERR030898.sickle.fq.gz


[ $? -ne 0 ] && echo `date` prostate ERR030898 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [prostate] [ERR030898] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030877

cd /home/wangq/data/bodymap2/process/prostate/ERR030877/

echo "* Start cufflinks [prostate] [ERR030877] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/prostate/ERR030877/cl_out \
    /home/wangq/data/bodymap2/process/prostate/ERR030877/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` prostate ERR030877 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [prostate] [ERR030877] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030898

cd /home/wangq/data/bodymap2/process/prostate/ERR030898/

echo "* Start cufflinks [prostate] [ERR030898] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/prostate/ERR030898/cl_out \
    /home/wangq/data/bodymap2/process/prostate/ERR030898/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` prostate ERR030898 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [prostate] [ERR030898] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


