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
# lane ERR030875

cd /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/

echo "* Start tophat [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (pair end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_1.sickle.fq.gz \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/trimmed/ERR030875_2.sickle.fq.gz


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# lane ERR030900

cd /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/

echo "* Start tophat [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log

# tophat (single end)
/home/wangq/bin/tophat -p 8 \
    -G /home/wangq/data/bodymap2/ref/human.37.gtf \
    -o /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/th_out \
    /home/wangq/data/bodymap2/ref/human.37 \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/trimmed/ERR030900.sickle.fq.gz


[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [tophat] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End tophat [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq/data/bodymap2/log/tophat.log


#----------------------------#
# cufflinks
#----------------------------#
# lane ERR030875

cd /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/

echo "* Start cufflinks [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/cl_out \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` white_blood_cells ERR030875 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [white_blood_cells] [ERR030875] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# lane ERR030900

cd /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/

echo "* Start cufflinks [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log

# cufflinks
/home/wangq/bin/cufflinks -p 8 \
    --no-update-check \
    -o /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/cl_out \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/th_out/accepted_hits.bam

[ $? -ne 0 ] && echo `date` white_blood_cells ERR030900 [cufflinks] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cufflinks [white_blood_cells] [ERR030900] `date`" | tee -a /home/wangq/data/bodymap2/log/cufflinks.log


