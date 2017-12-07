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
# lane ERR030864

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/

echo "* Start tophat [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/trimmed/ERR030864.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030865

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/

echo "* Start tophat [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/trimmed/ERR030865.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030857

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/

echo "* Start tophat [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/trimmed/ERR030857.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030858

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/

echo "* Start tophat [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/trimmed/ERR030858.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030856

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/

echo "* Start tophat [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/trimmed/ERR030856.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030864

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/

echo "* Start cufflinks [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030864/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030864 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_3] [ERR030864] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030865

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/

echo "* Start cufflinks [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030865/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030865 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_3] [ERR030865] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030857

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/

echo "* Start cufflinks [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030857/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030857 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_3] [ERR030857] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030858

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/

echo "* Start cufflinks [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030858/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030858 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_3] [ERR030858] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030856

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/

echo "* Start cufflinks [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_3/ERR030856/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_3 ERR030856 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_3] [ERR030856] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


