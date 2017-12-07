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
# lane ERR030888

cd /home/wangq/data/bodymap2/process/adipose/ERR030888/

echo "* Start tophat [adipose] [ERR030888] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/adipose/ERR030888/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/adipose/ERR030888/trimmed/ERR030888.sickle.fq.gz


[ $? -ne 0 ] && echo `date` adipose ERR030888 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [adipose] [ERR030888] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030880

cd /home/wangq/data/bodymap2/process/adipose/ERR030880/

echo "* Start tophat [adipose] [ERR030880] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/adipose/ERR030880/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/adipose/ERR030880/trimmed/ERR030880_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` adipose ERR030880 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [adipose] [ERR030880] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030888

cd /home/wangq/data/bodymap2/process/adipose/ERR030888/

echo "* Start cufflinks [adipose] [ERR030888] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/adipose/ERR030888/cl_out \
    /home/wangq/data/bodymap2/process/adipose/ERR030888/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` adipose ERR030888 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [adipose] [ERR030888] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030880

cd /home/wangq/data/bodymap2/process/adipose/ERR030880/

echo "* Start cufflinks [adipose] [ERR030880] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/adipose/ERR030880/cl_out \
    /home/wangq/data/bodymap2/process/adipose/ERR030880/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` adipose ERR030880 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [adipose] [ERR030880] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


