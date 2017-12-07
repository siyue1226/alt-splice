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

cd /home/wangq/data/bodymap2/process

#----------------------------#
# config file
#----------------------------#
if [ -e /home/wangq/data/bodymap2/process/assemblies.txt ]
then
    rm /home/wangq/data/bodymap2/process/assemblies.txt
fi

echo /home/wangq/data/bodymap2/process/adipose/ERR030888/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/adipose/ERR030880/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/adrenal/ERR030889/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/adrenal/ERR030881/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/brain/ERR030890/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/brain/ERR030882/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/breast/ERR030891/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/breast/ERR030883/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/colon/ERR030892/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/colon/ERR030884/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/heart/ERR030894/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/heart/ERR030886/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/kidney/ERR030885/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/kidney/ERR030893/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/liver/ERR030895/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/liver/ERR030887/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/lung/ERR030896/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/lung/ERR030879/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/lymph_node/ERR030897/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/lymph_node/ERR030878/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/ovary/ERR030874/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/ovary/ERR030901/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/prostate/ERR030877/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/prostate/ERR030898/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/testes/ERR030902/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/testes/ERR030873/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/thyroid/ERR030872/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/thyroid/ERR030903/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt
echo /home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/cl_out/transcripts.gtf >> /home/wangq/data/bodymap2/process/assemblies.txt

#----------------------------#
# cuffmerge
#----------------------------#
if [ -d /home/wangq/data/bodymap2/process/merged_asm ]
then
    rm -fr /home/wangq/data/bodymap2/process/merged_asm
fi

echo "* Start cuffmerge `date`" | tee -a /home/wangq/data/bodymap2/log/cuffmergediff.log

/home/wangq/bin/cuffmerge -p 14 \
    -g /home/wangq/data/bodymap2/ref/human.37.gtf \
    -s /home/wangq/data/bodymap2/ref/human.37.fa \
    -o /home/wangq/data/bodymap2/process/merged_asm \
    /home/wangq/data/bodymap2/process/assemblies.txt

[ $? -ne 0 ] && echo `date` ] [cuffmerge] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cuffmerge `date`" | tee -a /home/wangq/data/bodymap2/log/cuffmergediff.log


cd /home/wangq/data/bodymap2/process

#----------------------------#
# cuffdiff
#----------------------------#
echo "* Start cuffdiff `date`" | tee -a /home/wangq/data/bodymap2/log/cuffmergediff.log

/home/wangq/bin/cuffdiff -p 14 \
    --no-update-check \
    -u /home/wangq/data/bodymap2/process/merged_asm/merged.gtf \
    -b /home/wangq/data/bodymap2/ref/human.37.fa \
    -L adipose,adrenal,brain,breast,colon,heart,kidney,liver,lung,lymph_node,ovary,prostate,skeletal_muscle,testes,thyroid,white_blood_cells \
    /home/wangq/data/bodymap2/process/adipose/ERR030888/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/adipose/ERR030880/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/adrenal/ERR030889/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/adrenal/ERR030881/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/brain/ERR030890/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/brain/ERR030882/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/breast/ERR030891/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/breast/ERR030883/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/colon/ERR030892/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/colon/ERR030884/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/heart/ERR030894/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/heart/ERR030886/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/kidney/ERR030885/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/kidney/ERR030893/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/liver/ERR030895/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/liver/ERR030887/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/lung/ERR030896/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/lung/ERR030879/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/lymph_node/ERR030897/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/lymph_node/ERR030878/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/ovary/ERR030874/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/ovary/ERR030901/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/prostate/ERR030877/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/prostate/ERR030898/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/skeletal_muscle/ERR030876/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/skeletal_muscle/ERR030899/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/testes/ERR030902/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/testes/ERR030873/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/thyroid/ERR030872/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/thyroid/ERR030903/th_out/accepted_hits.bam \
    /home/wangq/data/bodymap2/process/white_blood_cells/ERR030875/th_out/accepted_hits.bam,/home/wangq/data/bodymap2/process/white_blood_cells/ERR030900/th_out/accepted_hits.bam \
    -o /home/wangq/data/bodymap2/process/diff_out

[ $? -ne 0 ] && echo `date` ] [cuffdiff] failed >> /home/wangq/data/bodymap2/fail.log && exit 255
echo "* End cuffdiff `date`" | tee -a /home/wangq/data/bodymap2/log/cuffmergediff.log

