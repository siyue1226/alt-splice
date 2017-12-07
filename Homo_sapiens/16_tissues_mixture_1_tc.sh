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
# lane ERR030868

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/

echo "* Start tophat [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/trimmed/ERR030868.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030869

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/

echo "* Start tophat [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/trimmed/ERR030869.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030871

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/

echo "* Start tophat [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/trimmed/ERR030871.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030863

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/

echo "* Start tophat [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/trimmed/ERR030863.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030862

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/

echo "* Start tophat [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/trimmed/ERR030862.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030870

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/

echo "* Start tophat [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/trimmed/ERR030870.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030868

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030868/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030868 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030868] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030869

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030869/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030869 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030869] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030871

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030871/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030871 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030871] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030863

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030863/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030863 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030863] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030862

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030862/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030862 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030862] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030870

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/

echo "* Start cufflinks [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_1/ERR030870/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_1 ERR030870 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_1] [ERR030870] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


