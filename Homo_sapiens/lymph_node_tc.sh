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
# lane ERR030897

cd /home/wangq/data/bodymap2/process/lymph_node/ERR030897/

echo "* Start tophat [lymph_node] [ERR030897] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/lymph_node/ERR030897/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030897/trimmed/ERR030897.sickle.fq.gz


[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [lymph_node] [ERR030897] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030878

cd /home/wangq/data/bodymap2/process/lymph_node/ERR030878/

echo "* Start tophat [lymph_node] [ERR030878] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/lymph_node/ERR030878/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030878/trimmed/ERR030878_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [lymph_node] [ERR030878] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030897

cd /home/wangq/data/bodymap2/process/lymph_node/ERR030897/

echo "* Start cufflinks [lymph_node] [ERR030897] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/lymph_node/ERR030897/cl_out \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030897/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` lymph_node ERR030897 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [lymph_node] [ERR030897] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030878

cd /home/wangq/data/bodymap2/process/lymph_node/ERR030878/

echo "* Start cufflinks [lymph_node] [ERR030878] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/lymph_node/ERR030878/cl_out \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030878/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` lymph_node ERR030878 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [lymph_node] [ERR030878] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


