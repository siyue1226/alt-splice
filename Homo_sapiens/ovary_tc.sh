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
# lane ERR030874

cd /home/wangq/data/bodymap2/process/ovary/ERR030874/

echo "* Start tophat [ovary] [ERR030874] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/ovary/ERR030874/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/ovary/ERR030874/trimmed/ERR030874_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` ovary ERR030874 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [ovary] [ERR030874] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030901

cd /home/wangq/data/bodymap2/process/ovary/ERR030901/

echo "* Start tophat [ovary] [ERR030901] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/ovary/ERR030901/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/ovary/ERR030901/trimmed/ERR030901.sickle.fq.gz


[ $? -ne 0 ] && echo `date` ovary ERR030901 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [ovary] [ERR030901] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030874

cd /home/wangq/data/bodymap2/process/ovary/ERR030874/

echo "* Start cufflinks [ovary] [ERR030874] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/ovary/ERR030874/cl_out \
    /home/wangq/data/bodymap2/process/ovary/ERR030874/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` ovary ERR030874 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [ovary] [ERR030874] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030901

cd /home/wangq/data/bodymap2/process/ovary/ERR030901/

echo "* Start cufflinks [ovary] [ERR030901] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/ovary/ERR030901/cl_out \
    /home/wangq/data/bodymap2/process/ovary/ERR030901/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` ovary ERR030901 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [ovary] [ERR030901] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


