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
# lane ERR030876

cd /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/

echo "* Start tophat [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/trimmed/ERR030876_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030899

cd /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/

echo "* Start tophat [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/trimmed/ERR030899.sickle.fq.gz


[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030876

cd /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/

echo "* Start cufflinks [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/cl_out \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030876 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [skeletal_muscle] [ERR030876] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030899

cd /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/

echo "* Start cufflinks [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/cl_out \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` skeletal_muscle ERR030899 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [skeletal_muscle] [ERR030899] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


