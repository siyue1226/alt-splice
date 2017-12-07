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
# lane ERR030859

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/

echo "* Start tophat [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/trimmed/ERR030859.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030867

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/

echo "* Start tophat [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/trimmed/ERR030867.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030861

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/

echo "* Start tophat [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/trimmed/ERR030861.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030866

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/

echo "* Start tophat [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/trimmed/ERR030866.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030860

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/

echo "* Start tophat [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/trimmed/ERR030860.sickle.fq.gz


[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030859

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/

echo "* Start cufflinks [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030859/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030859 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_2] [ERR030859] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030867

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/

echo "* Start cufflinks [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030867/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030867 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_2] [ERR030867] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030861

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/

echo "* Start cufflinks [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030861/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030861 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_2] [ERR030861] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030866

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/

echo "* Start cufflinks [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030866/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030866 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_2] [ERR030866] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030860

cd /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/

echo "* Start cufflinks [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/cl_out \
    /home/wangq/data/bodymap2/process/16_tissues_mixture_2/ERR030860/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` 16_tissues_mixture_2 ERR030860 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [16_tissues_mixture_2] [ERR030860] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


