#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/rice_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa /home/wangq/data/rna-seq/rice_trans/ref/rice.29

cd /home/wangq/data/rna-seq/rice_trans/process

#----------------------------#
# config file
#----------------------------#
if [ -e /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt ]
then
    rm /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
fi

echo /home/wangq/data/rna-seq/rice_trans/process/SRR352184/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352187/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352189/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352190/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352192/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352194/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352204/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352206/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352207/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352209/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/rice_trans/process/SRR352211/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt



#----------------------------#
# cuffmerge
#----------------------------#
if [ -d /home/wangq/data/rna-seq/rice_trans/process/merged_asm ]
then
    rm -fr /home/wangq/data/rna-seq/rice_trans/process/merged_asm
fi

echo "* Start cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cuffmergediff.log

/home/wangq/bin/cuffmerge -p 14 \
    -g /home/wangq/data/rna-seq/rice_trans/ref/rice.29.gtf \
    -s /home/wangq/data/rna-seq/rice_trans/ref/rice.29.fa \
    -o /home/wangq/data/rna-seq/rice_trans/process/merged_asm \
    /home/wangq/data/rna-seq/rice_trans/process/assemblies.txt

[ $? -ne 0 ] && echo `date` ] [cuffmerge] failed >> /home/wangq/data/rna-seq/rice_trans/fail.log && exit 255
echo "* End cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/rice_trans/log/cuffmergediff.log


