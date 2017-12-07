#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/mouse_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65

cd /home/wangq/data/rna-seq/mouse_trans/process

#----------------------------#
# config file
#----------------------------#
if [ -e /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt ]
then
    rm /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
fi

echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453116/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453117/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453118/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453119/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453120/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Adrenal/SRR453121/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453166/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453167/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453168/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453169/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453170/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Colon/SRR453171/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453109/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453110/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453111/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453112/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453113/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453114/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Duodenum/SRR453115/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453126/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453127/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453128/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/GenitalFatPad/SRR453129/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453172/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453173/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453174/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Heart/SRR453175/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453144/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453145/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453146/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453147/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453148/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Kidney/SRR453149/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453122/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453123/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453124/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/LgIntestine/SRR453125/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453150/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453151/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453152/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453153/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453154/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Liver/SRR453155/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453156/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453157/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453158/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Lung/SRR453159/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453087/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453088/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453089/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453090/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453091/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/MammaryGland/SRR453092/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453077/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453078/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453079/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453080/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453081/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453082/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453083/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453084/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453085/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Ovary/SRR453086/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453099/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453100/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453101/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453102/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453103/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453104/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453105/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453106/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453107/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SmIntestine/SRR453108/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453160/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453161/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453162/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453163/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453164/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Spleen/SRR453165/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453093/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453094/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453095/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453096/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453097/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Stomach/SRR453098/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453130/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453131/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453132/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/SubcFatPad/SRR453133/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453140/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453141/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453142/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Testis/SRR453143/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453134/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453135/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453136/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453137/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453138/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/mouse_trans/process/Thymus/SRR453139/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt

#----------------------------#
# cuffmerge
#----------------------------#
if [ -d /home/wangq/data/rna-seq/mouse_trans/process/merged_asm ]
then
    rm -fr /home/wangq/data/rna-seq/mouse_trans/process/merged_asm
fi

echo "* Start cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffmergediff.log

/home/wangq/bin/cuffmerge -p 14 \
    -g /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.gtf \
    -s /home/wangq/data/rna-seq/mouse_trans/ref/mouse.65.fa \
    -o /home/wangq/data/rna-seq/mouse_trans/process/merged_asm \
    /home/wangq/data/rna-seq/mouse_trans/process/assemblies.txt

[ $? -ne 0 ] && echo `date` ] [cuffmerge] failed >> /home/wangq/data/rna-seq/mouse_trans/fail.log && exit 255
echo "* End cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/mouse_trans/log/cuffmergediff.log


