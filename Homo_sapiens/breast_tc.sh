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
# lane ERR030891

cd /home/wangq/data/bodymap2/process/breast/ERR030891/

echo "* Start tophat [breast] [ERR030891] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/breast/ERR030891/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/breast/ERR030891/trimmed/ERR030891.sickle.fq.gz


[ $? -ne 0 ] && echo `date` breast ERR030891 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [breast] [ERR030891] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030883

cd /home/wangq/data/bodymap2/process/breast/ERR030883/

echo "* Start tophat [breast] [ERR030883] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/breast/ERR030883/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/breast/ERR030883/trimmed/ERR030883_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/breast/ERR030883/trimmed/ERR030883_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` breast ERR030883 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [breast] [ERR030883] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030891

cd /home/wangq/data/bodymap2/process/breast/ERR030891/

echo "* Start cufflinks [breast] [ERR030891] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/breast/ERR030891/cl_out \
    /home/wangq/data/bodymap2/process/breast/ERR030891/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` breast ERR030891 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [breast] [ERR030891] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030883

cd /home/wangq/data/bodymap2/process/breast/ERR030883/

echo "* Start cufflinks [breast] [ERR030883] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/breast/ERR030883/cl_out \
    /home/wangq/data/bodymap2/process/breast/ERR030883/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` breast ERR030883 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [breast] [ERR030883] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


