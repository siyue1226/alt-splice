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
# lane ERR030872

cd /home/wangq/data/bodymap2/process/thyroid/ERR030872/

echo "* Start tophat [thyroid] [ERR030872] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/thyroid/ERR030872/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/thyroid/ERR030872/trimmed/ERR030872_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` thyroid ERR030872 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [thyroid] [ERR030872] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030903

cd /home/wangq/data/bodymap2/process/thyroid/ERR030903/

echo "* Start tophat [thyroid] [ERR030903] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/thyroid/ERR030903/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/thyroid/ERR030903/trimmed/ERR030903.sickle.fq.gz


[ $? -ne 0 ] && echo `date` thyroid ERR030903 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [thyroid] [ERR030903] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030872

cd /home/wangq/data/bodymap2/process/thyroid/ERR030872/

echo "* Start cufflinks [thyroid] [ERR030872] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/thyroid/ERR030872/cl_out \
    /home/wangq/data/bodymap2/process/thyroid/ERR030872/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` thyroid ERR030872 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [thyroid] [ERR030872] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030903

cd /home/wangq/data/bodymap2/process/thyroid/ERR030903/

echo "* Start cufflinks [thyroid] [ERR030903] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/thyroid/ERR030903/cl_out \
    /home/wangq/data/bodymap2/process/thyroid/ERR030903/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` thyroid ERR030903 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [thyroid] [ERR030903] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


