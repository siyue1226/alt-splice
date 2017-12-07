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
# lane ERR030885

cd /home/wangq/data/bodymap2/process/kidney/ERR030885/

echo "* Start tophat [kidney] [ERR030885] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/kidney/ERR030885/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/kidney/ERR030885/trimmed/ERR030885_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` kidney ERR030885 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [kidney] [ERR030885] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030893

cd /home/wangq/data/bodymap2/process/kidney/ERR030893/

echo "* Start tophat [kidney] [ERR030893] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/kidney/ERR030893/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/kidney/ERR030893/trimmed/ERR030893.sickle.fq.gz


[ $? -ne 0 ] && echo `date` kidney ERR030893 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [kidney] [ERR030893] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030885

cd /home/wangq/data/bodymap2/process/kidney/ERR030885/

echo "* Start cufflinks [kidney] [ERR030885] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/kidney/ERR030885/cl_out \
    /home/wangq/data/bodymap2/process/kidney/ERR030885/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` kidney ERR030885 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [kidney] [ERR030885] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030893

cd /home/wangq/data/bodymap2/process/kidney/ERR030893/

echo "* Start cufflinks [kidney] [ERR030893] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/kidney/ERR030893/cl_out \
    /home/wangq/data/bodymap2/process/kidney/ERR030893/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` kidney ERR030893 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [kidney] [ERR030893] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


