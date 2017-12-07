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
# cuffdiff_cxb
#----------------------------#
echo "* Start cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffmergediff.log

/home/wangq/bin/cuffdiff -p 14 \
    --no-update-check -u -M /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.mask.gtf \
    -b /home/wangq/data/rna-seq/dmel_trans/ref/dmel.12_65.fa \
    /home/wangq/data/rna-seq/dmel_trans/process/merged_asm/merged.gtf \
    -L WPP_2d_CNS,WPP_2d_fat_body,WPP_fat_body,WPP_salivary_glands,mated_female_eclosion_1d_heads,mated_female_eclosion_20d_heads,mated_female_eclosion_4d_heads,mated_female_eclosion_4d_ovaries,mated_male_eclosion_1d_heads,mated_male_eclosion_20d_heads,mated_male_eclosion_4d_accessory_glands,mated_male_eclosion_4d_heads,mated_male_eclosion_4d_testes,mixed_males_females_eclosion_1d_carcass,mixed_males_females_eclosion_1d_digestive_system,mixed_males_females_eclosion_20d_carcass,mixed_males_females_eclosion_20d_digestive_system,mixed_males_females_eclosion_4d_carcass,mixed_males_females_eclosion_4d_digestive_system,third_instar_larvae_wandering_stage_CNS,third_instar_larvae_wandering_stage_carcass,third_instar_larvae_wandering_stage_digestive_system,third_instar_larvae_wandering_stage_fat_body,third_instar_larvae_wandering_stage_imaginal_discs,third_instar_larvae_wandering_stage_salivary_glands,virgin_female_eclosion_1d_heads,virgin_female_eclosion_20d_heads,virgin_female_eclosion_4d_heads,virgin_female_eclosion_4d_ovaries \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR070412/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_CNS/SRR100271/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_2d_fat_body/SRR070429/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070411/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/WPP_fat_body/SRR070428/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR070427/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/WPP_salivary_glands/SRR100270/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070434/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR070435/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_1d_heads/SRR100279/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR070420/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR100274/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR116383/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_20d_heads/SRR111882/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070414/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_heads/SRR070415/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR070431/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100277/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_female_eclosion_4d_ovaries/SRR100283/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070432/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR070433/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_1d_heads/SRR100280/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070421/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_20d_heads/SRR070424/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR070397/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR182358/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_accessory_glands/SRR350959/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070400/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_heads/SRR070416/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070422/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR070423/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR100276/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350960/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mated_male_eclosion_4d_testes/SRR350961/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070399/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_carcass/SRR070395/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070398/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_1d_digestive_system/SRR070394/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070404/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_carcass/SRR070391/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR070403/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_20d_digestive_system/SRR111883/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_carcass/SRR070387/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_carcass/SRR070402/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR070401/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111878/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/mixed_males_females_eclosion_4d_digestive_system/SRR111879/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070409/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_CNS/SRR070410/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR070426/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_carcass/SRR100269/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR070408/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_digestive_system/SRR100268/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070405/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_fat_body/SRR070406/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070392/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR070393/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111884/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR111885/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350962/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_imaginal_discs/SRR350963/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070407/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/third_instar_larvae_wandering_stage_salivary_glands/SRR070425/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070436/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR070437/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_1d_heads/SRR100281/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070388/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR070419/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_20d_heads/SRR100275/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR070430/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100278/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_heads/SRR100282/cq_out/abundances.cxb \
    /home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070396/cq_out/abundances.cxb,/home/wangq/data/rna-seq/dmel_trans/process/virgin_female_eclosion_4d_ovaries/SRR070417/cq_out/abundances.cxb \
    -o /home/wangq/data/rna-seq/dmel_trans/process/diff_out_cxb

[ $? -ne 0 ] && echo `date` ] [cuffdiff_cxb] failed >> /home/wangq/data/rna-seq/dmel_trans/fail.log && exit 255
echo "* End cuffdiff_cxb `date`" | tee -a /home/wangq/data/rna-seq/dmel_trans/log/cuffmergediff.log

