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
# lane ERR030896

cd /home/wangq/data/bodymap2/process/lung/ERR030896/

echo "* Start tophat [lung] [ERR030896] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/lung/ERR030896/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/lung/ERR030896/trimmed/ERR030896.sickle.fq.gz


[ $? -ne 0 ] && echo `date` lung ERR030896 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [lung] [ERR030896] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030879

cd /home/wangq/data/bodymap2/process/lung/ERR030879/

echo "* Start tophat [lung] [ERR030879] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/lung/ERR030879/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/lung/ERR030879/trimmed/ERR030879_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/lung/ERR030879/trimmed/ERR030879_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` lung ERR030879 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [lung] [ERR030879] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030896

cd /home/wangq/data/bodymap2/process/lung/ERR030896/

echo "* Start cufflinks [lung] [ERR030896] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/lung/ERR030896/cl_out \
    /home/wangq/data/bodymap2/process/lung/ERR030896/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` lung ERR030896 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [lung] [ERR030896] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030879

cd /home/wangq/data/bodymap2/process/lung/ERR030879/

echo "* Start cufflinks [lung] [ERR030879] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/lung/ERR030879/cl_out \
    /home/wangq/data/bodymap2/process/lung/ERR030879/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` lung ERR030879 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [lung] [ERR030879] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


