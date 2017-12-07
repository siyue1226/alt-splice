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
# lane ERR030895

cd /home/wangq/data/bodymap2/process/liver/ERR030895/

echo "* Start tophat [liver] [ERR030895] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/liver/ERR030895/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/liver/ERR030895/trimmed/ERR030895.sickle.fq.gz


[ $? -ne 0 ] && echo `date` liver ERR030895 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [liver] [ERR030895] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030887

cd /home/wangq/data/bodymap2/process/liver/ERR030887/

echo "* Start tophat [liver] [ERR030887] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/liver/ERR030887/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/liver/ERR030887/trimmed/ERR030887_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/liver/ERR030887/trimmed/ERR030887_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` liver ERR030887 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [liver] [ERR030887] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030895

cd /home/wangq/data/bodymap2/process/liver/ERR030895/

echo "* Start cufflinks [liver] [ERR030895] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/liver/ERR030895/cl_out \
    /home/wangq/data/bodymap2/process/liver/ERR030895/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` liver ERR030895 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [liver] [ERR030895] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030887

cd /home/wangq/data/bodymap2/process/liver/ERR030887/

echo "* Start cufflinks [liver] [ERR030887] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/liver/ERR030887/cl_out \
    /home/wangq/data/bodymap2/process/liver/ERR030887/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` liver ERR030887 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [liver] [ERR030887] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


