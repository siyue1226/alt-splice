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
# lane ERR030902

cd /home/wangq/data/bodymap2/process/testes/ERR030902/

echo "* Start tophat [testes] [ERR030902] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/testes/ERR030902/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/testes/ERR030902/trimmed/ERR030902.sickle.fq.gz


[ $? -ne 0 ] && echo `date` testes ERR030902 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [testes] [ERR030902] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030873

cd /home/wangq/data/bodymap2/process/testes/ERR030873/

echo "* Start tophat [testes] [ERR030873] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/testes/ERR030873/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/testes/ERR030873/trimmed/ERR030873_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/testes/ERR030873/trimmed/ERR030873_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` testes ERR030873 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [testes] [ERR030873] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030902

cd /home/wangq/data/bodymap2/process/testes/ERR030902/

echo "* Start cufflinks [testes] [ERR030902] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/testes/ERR030902/cl_out \
    /home/wangq/data/bodymap2/process/testes/ERR030902/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` testes ERR030902 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [testes] [ERR030902] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030873

cd /home/wangq/data/bodymap2/process/testes/ERR030873/

echo "* Start cufflinks [testes] [ERR030873] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/testes/ERR030873/cl_out \
    /home/wangq/data/bodymap2/process/testes/ERR030873/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` testes ERR030873 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [testes] [ERR030873] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


