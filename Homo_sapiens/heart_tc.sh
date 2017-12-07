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
# lane ERR030894

cd /home/wangq/data/bodymap2/process/heart/ERR030894/

echo "* Start tophat [heart] [ERR030894] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/heart/ERR030894/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/heart/ERR030894/trimmed/ERR030894.sickle.fq.gz


[ $? -ne 0 ] && echo `date` heart ERR030894 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [heart] [ERR030894] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030886

cd /home/wangq/data/bodymap2/process/heart/ERR030886/

echo "* Start tophat [heart] [ERR030886] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/heart/ERR030886/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/heart/ERR030886/trimmed/ERR030886_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/heart/ERR030886/trimmed/ERR030886_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` heart ERR030886 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [heart] [ERR030886] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030894

cd /home/wangq/data/bodymap2/process/heart/ERR030894/

echo "* Start cufflinks [heart] [ERR030894] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/heart/ERR030894/cl_out \
    /home/wangq/data/bodymap2/process/heart/ERR030894/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` heart ERR030894 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [heart] [ERR030894] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030886

cd /home/wangq/data/bodymap2/process/heart/ERR030886/

echo "* Start cufflinks [heart] [ERR030886] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/heart/ERR030886/cl_out \
    /home/wangq/data/bodymap2/process/heart/ERR030886/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` heart ERR030886 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [heart] [ERR030886] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


