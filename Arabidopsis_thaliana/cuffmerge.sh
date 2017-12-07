#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/Ath_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29

cd /home/wangq/data/rna-seq/Ath_trans/process

#----------------------------#
# config file
#----------------------------#
if [ -e /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt ]
then
    rm /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt
fi

echo /home/wangq/data/rna-seq/Ath_trans/process/SRR6179907/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/Ath_trans/process/SRR6179906/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/Ath_trans/process/SRR6179904/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/Ath_trans/process/SRR6179905/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/Ath_trans/process/SRR6179908/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt


#----------------------------#
# cuffmerge
#----------------------------#
if [ -d /home/wangq/data/rna-seq/Ath_trans/process/merged_asm ]
then
    rm -fr /home/wangq/data/rna-seq/Ath_trans/process/merged_asm
fi

echo "* Start cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffmergediff.log

/home/wangq/bin/cuffmerge -p 14 \
    -g /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.gtf \
    -s /home/wangq/data/rna-seq/Ath_trans/ref/Ath.29.fa \
    -o /home/wangq/data/rna-seq/Ath_trans/process/merged_asm \
    /home/wangq/data/rna-seq/Ath_trans/process/assemblies.txt

[ $? -ne 0 ] && echo `date` ] [cuffmerge] failed >> /home/wangq/data/rna-seq/Ath_trans/fail.log && exit 255
echo "* End cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/Ath_trans/log/cuffmergediff.log


