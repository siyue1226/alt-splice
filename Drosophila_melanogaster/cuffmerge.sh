#!/bin/bash
start_time=`date +%s`

### bash hints
### 2>&1                        redirect stderr to stdout
### | tee -a log.log            screen outputs also append to log file
### ; ( exit ${PIPESTATUS} )    correct program exitting status
### Only run parallel when you're sure that there are no errors.

cd /home/wangq/data/rna-seq/dmel_trans

### bowtie index
# bowtie2-build /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65

cd /home/wangq/data/rna-seq/dmel_trans/process

#----------------------------#
# config file
#----------------------------#
if [ -e /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt ]
then
    rm /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
fi

echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_carcass/SRR070387/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_carcass/SRR070402/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070405/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070406/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt
echo /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/cl_out/transcripts.gtf >> /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt

#----------------------------#
# cuffmerge
#----------------------------#
if [ -d /home/wangq/data/rna-seq/dmel_trans/process/merged_asm ]
then
    rm -fr /home/wangq/data/rna-seq/dmel_trans/process/merged_asm
fi

echo "* Start cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffmergediff.log

/home/wangq/bin/cuffmerge -p 14 \
    -g /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.gtf \
    -s /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    -o /home/wangq/data/rna-seq/dmel_trans/process/merged_asm \
    /home/wangq/data/rna-seq/dmel_trans/process/assemblies.txt

[ $? -ne 0 ] && echo `date` ] [cuffmerge] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffmerge `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffmergediff.log


