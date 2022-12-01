#==================
#1. load in all files
#==================
alzheimers_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90027158_buildGRCh38.tsv', sep = '\t', header = TRUE)
body_mass_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/BMI_European.fmt.gzip', sep = '\t', header = TRUE)
chronic_inflammation_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90029070_buildGRCh37.tsv', sep = '\t', header = TRUE)
class_1_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS1_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
class_2_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS2_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
class_3_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS3_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
coronary_artery_disease_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/29212778-GCST005195-EFO_0000378-build37.f.tsv', sep = '\t', header = TRUE)
extreme_body_mass_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_EXTREME_BMI_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
extreme_waist_hip_ratio_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_EXTREME_WHR_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
extrinsic_epigenetic_age_acceleration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/EEAA_metaPlusGS.csv', sep = ',', header = TRUE)
fasting_glucose_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt', sep = '\t', header = TRUE)
fasting_insulin_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt', sep = '\t', header = TRUE)
fasting_proinsulin_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_ln_fastingProinsulin.txt', sep = '\t', header = TRUE)
frailty_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/Atkins_2020_Frailty_Index_GWAS_UKB_Europeans_60_to_70_years.txt', sep = '\t', header = TRUE)
hba1c_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_HbA1C.txt', sep = ' ', header = TRUE)
healthspan_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/healthspan_summary.csv', sep = ',', header = TRUE)
high_density_lipoprotien_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE)
homab_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_HOMAB_MainEffect.txt', sep = ' ', header = TRUE)
homair_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_HOMAIR_MainEffect.txt', sep = ' ', header = TRUE)
income_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/HillWD_31844048_household_Income.txt', sep = ' ', header = TRUE)
intelligence_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/29942086-GCST006250-EFO_0004337-Build37.f.tsv', sep = '\t', header = TRUE)
intrinsic_epigenetic_age_acceleration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/IEAA_metaPlusGS.csv', sep = ',', header = TRUE)
lacunar_stroke_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90014122_buildGRCh37.tsv', sep = '\t', header = TRUE)
low_density_lipoprotien_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE)
overweight_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OVERWEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE)
percieved_age_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/perceived_age_summary_results.txt', sep = ' ', header = TRUE)
serum_urate_concentration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/23263486-GCST001791-EFO_0004531-build36.f.tsv', sep = '\t', header = TRUE)
telomere_length_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/UKB_telomere_gwas_summarystats.tsv', sep = '\t', header = TRUE)
triglyceride_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE)
two_hour_glucose_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_2hrGlucose_AdjustedForBMI.txt', sep = '\t', header = TRUE)
type_2_diabetes_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/DIAGRAMv3.2012DEC17.txt', sep = '\t', header = TRUE)
waist_hip_ratio_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt', sep = ' ', header = TRUE)
atherosclerosis_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90044005_buildGRCh37.tsv', sep = '\t', header = TRUE)
leucine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132732_buildGRCh37.txt', sep = ' ', header = TRUE)
isoleucine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132730_buildGRCh37.txt', sep = ' ', header = TRUE)
valine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132752_buildGRCh37.txt', sep = ' ', header = TRUE)
free_cholesterol_in_very_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132680_buildGRCh37.txt', sep = ' ', header = TRUE)
free_cholesterol_in_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132677_buildGRCh37.txt', sep = ' ', header = TRUE)
phospholipid_in_very_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132642_buildGRCh37.txt', sep = ' ', header = TRUE)
phospholipid_in_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132640_buildGRCh37.txt', sep = ' ', header = TRUE)
triglyceride_in_very_large_vldl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132668_buildGRCh37.txt', sep = ' ', header = TRUE)
triglyceride_in_large_vldl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132665_buildGRCh37.txt', sep = ' ', header = TRUE)












# alzheimers_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90027158_buildGRCh38.tsv', sep = '\t', header = TRUE, nrows=10)
# body_mass_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/BMI_European.fmt.gzip', sep = '\t', header = TRUE, nrows=10)
# chronic_inflammation_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90029070_buildGRCh37.tsv', sep = '\t', header = TRUE, nrows=10)
# class_1_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS1_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# class_2_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS2_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# class_3_obesity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OBESITY_CLASS3_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# coronary_artery_disease_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/29212778-GCST005195-EFO_0000378-build37.f.tsv', sep = '\t', header = TRUE, nrows=10)
# extreme_body_mass_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_EXTREME_BMI_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# extreme_waist_hip_ratio_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_EXTREME_WHR_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# extrinsic_epigenetic_age_acceleration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/EEAA_metaPlusGS.csv', sep = ',', header = TRUE, nrows=10)
# fasting_glucose_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt', sep = '\t', header = TRUE, nrows=10)
# fasting_insulin_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt', sep = '\t', header = TRUE, nrows=10)
# fasting_proinsulin_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_ln_fastingProinsulin.txt', sep = '\t', header = TRUE, nrows=10)
# frailty_index_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/Atkins_2020_Frailty_Index_GWAS_UKB_Europeans_60_to_70_years.txt', sep = '\t', header = TRUE, nrows=10)
# hba1c_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_HbA1C.txt', sep = ' ', header = TRUE, nrows=10)
# healthspan_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/healthspan_summary.csv', sep = ',', header = TRUE, nrows=10)
# high_density_lipoprotien_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE, nrows=10)
# homab_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_HOMAB_MainEffect.txt', sep = ' ', header = TRUE, nrows=10)
# homair_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_Manning_et_al_HOMAIR_MainEffect.txt', sep = ' ', header = TRUE, nrows=10)
# income_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/HillWD_31844048_household_Income.txt', sep = ' ', header = TRUE, nrows=10)
# intelligence_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/29942086-GCST006250-EFO_0004337-Build37.f.tsv', sep = '\t', header = TRUE, nrows=10)
# intrinsic_epigenetic_age_acceleration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/IEAA_metaPlusGS.csv', sep = ',', header = TRUE, nrows=10)
# lacunar_stroke_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90014122_buildGRCh37.tsv', sep = '\t', header = TRUE, nrows=10)
# # longevity_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/Results_99th_percentile.txt', sep = '\t', header = TRUE, nrows=10)
# low_density_lipoprotien_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE, nrows=10)
# overweight_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GIANT_OVERWEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt', sep = ' ', header = TRUE, nrows=10)
# percieved_age_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/perceived_age_summary_results.txt', sep = ' ', header = TRUE, nrows=10)
# serum_urate_concentration_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/23263486-GCST001791-EFO_0004531-build36.f.tsv', sep = '\t', header = TRUE, nrows=10)
# telomere_length_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/UKB_telomere_gwas_summarystats.tsv', sep = '\t', header = TRUE, nrows=10)
# triglyceride_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results', sep = '\t', header = TRUE, nrows=10)
# two_hour_glucose_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/MAGIC_2hrGlucose_AdjustedForBMI.txt', sep = '\t', header = TRUE, nrows=10)
# type_2_diabetes_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/DIAGRAMv3.2012DEC17.txt', sep = '\t', header = TRUE, nrows=10)
# waist_hip_ratio_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt', sep = ' ', header = TRUE, nrows=10)
# atherosclerosis_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90044005_buildGRCh37.tsv', sep = '\t', header = TRUE, nrows=10)
# leucine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132732_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# isoleucine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132730_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# valine_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132752_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# free_cholesterol_in_very_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132680_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# free_cholesterol_in_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132677_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# phospholipid_in_very_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132642_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# phospholipid_in_large_hdl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132640_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# triglyceride_in_very_large_vldl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132668_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)
# triglyceride_in_large_vldl_master = read.table(file = '~/scratch/principal_component_GWAS/all_traits_raw/GCST90132665_buildGRCh37.txt', sep = ' ', header = TRUE, nrows=10)


print("1. load in all files")
#==================
#2. calculate and filter all columns
#==================
#alzheimers: n
alzheimers_master$n = alzheimers_master$n_cases + alzheimers_master$n_controls
alzheimers = alzheimers_master[,c("variant_id","other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

#body_mass_index: n
body_mass_index_master$n = rep(718734, nrow(body_mass_index_master))
body_mass_index = body_mass_index_master[,c("SNPNAME", "REF", "ALT", "n", "beta", "Pvalue", "se")]

#chronic_inflammation: n
chronic_inflammation_master$n = rep(1002898, nrow(chronic_inflammation_master))
chronic_inflammation = chronic_inflammation_master[,c("variant_id", "other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

#class_1_obesity: n
class_1_obesity_master$n = class_1_obesity_master$N_cases + class_1_obesity_master$N_controls
class_1_obesity = class_1_obesity_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#class_2_obesity: n

class_2_obesity_master$n = class_2_obesity_master$N_cases + class_2_obesity_master$N_controls
class_2_obesity = class_2_obesity_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#class_3_obesity: n
class_3_obesity_master$n = class_3_obesity_master$N_cases + class_3_obesity_master$N_controls
class_3_obesity = class_3_obesity_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#coronary_artery_disease: n/a
coronary_artery_disease = coronary_artery_disease_master[,c("variant_id", "other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

#extreme_body_mass_index: n
extreme_body_mass_index_master$n = extreme_body_mass_index_master$N_cases + extreme_body_mass_index_master$N_controls
extreme_body_mass_index = extreme_body_mass_index_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#extreme_waist_hip_ratio: n
extreme_waist_hip_ratio_master$n = extreme_waist_hip_ratio_master$N_cases + extreme_waist_hip_ratio_master$N_controls
extreme_waist_hip_ratio = extreme_waist_hip_ratio_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#extrinsic_epigenetic_age_acceleration: n
extrinsic_epigenetic_age_acceleration_master$n = rep(13493, nrow(extrinsic_epigenetic_age_acceleration_master))
extrinsic_epigenetic_age_acceleration = extrinsic_epigenetic_age_acceleration_master[,c("SNP", "metaA1", "metaA2", "n", "meta_b", "metaPval", "meta_se")]

#fasting_glucose: n
fasting_glucose_master$n = rep(58074, nrow(fasting_glucose_master))
fasting_glucose = fasting_glucose_master[,c("Snp", "other_allele", "effect_allele", "n", "BMIadjMainEffects", "BMIadjMainP", "BMIadjMainSE")]

#fasting_insulin: n
fasting_insulin_master$n = rep(51750, nrow(fasting_insulin_master))
fasting_insulin = fasting_insulin_master[,c("Snp", "other_allele", "effect_allele", "n", "BMIadjMainEffects", "BMIadjMainP", "BMIadjMainSE")]

#fasting_proinsulin: n
fasting_proinsulin_master$n = rep(10701, nrow(fasting_proinsulin_master))
fasting_proinsulin = fasting_proinsulin_master[,c("snp", "other_allele", "effect_allele", "n", "effect", "pvalue", "stderr")]

#frailty_index: n
frailty_index_master$n = rep(175226, nrow(frailty_index_master))
frailty_index = frailty_index_master[,c("SNP", "ALLELE0", "ALLELE1", "n", "BETA", "P_BOLT_LMM", "SE")]

#hba1c: n
hba1c_master$n = rep(46368, nrow(hba1c_master))
hba1c = hba1c_master[,c("snp", "other_allele", "effect_allele", "n", "effect", "pvalue", "stderr")]

#healthspan: n, p_value
healthspan_master$n = rep(300447, nrow(healthspan_master))
healthspan_master$p_value = 10^(-1*(healthspan_master$X.log10.p.value.))
healthspan = healthspan_master[,c("SNPID", "RA", "EA", "n", "beta", "p_value", "se")]

#high_density_lipoprotien: n/a
high_density_lipoprotien = high_density_lipoprotien_master[,c("rsID", "REF", "ALT", "N", "EFFECT_SIZE", "pvalue", "SE")]

#homab: n
homab_master$n = rep(58074, nrow(homab_master))
homab = homab_master[,c("Snp", "other_allele", "effect_allele", "n", "BMIadjMainEffects", "BMIadjMainP", "BMIadjMainSE")]

#homab: n
homair_master$n = rep(51750, nrow(homair_master))
homair = homair_master[,c("Snp", "other_allele", "effect_allele", "n", "BMIadjMainEffects", "BMIadjMainP", "BMIadjMainSE")]

#income_master: n
income_master$n = rep(286301, nrow(income_master))
income = income_master[,c("SNP", "Non_effect_Allele", "Effect_Allele", "n", "Beta", "P", "Standard_Error_of_Beta")]

#intelligence: n
intelligence_master$n = rep(269867, nrow(intelligence_master))
intelligence = intelligence_master[,c("variant_id", "other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

#intrinsic_epigenetic_age_acceleration: n
intrinsic_epigenetic_age_acceleration_master$n = rep(13493, nrow(intrinsic_epigenetic_age_acceleration_master))
intrinsic_epigenetic_age_acceleration = intrinsic_epigenetic_age_acceleration_master[,c("SNP", "metaA1", "metaA2", "n", "meta_b", "metaPval", "meta_se")]

#lacunar_stroke: n
lacunar_stroke_master$n = rep(254959, nrow(lacunar_stroke_master))
lacunar_stroke = lacunar_stroke_master[,c("variant_id", "other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

# #longevity: rsid
# library(biomaRt)
# snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP", 
#                        host="grch37.ensembl.org", 
#                        dataset="hsapiens_snp")
# positions_longevity = paste0(longevity_master$Chr, ":", longevity_master$Position, ":", longevity_master$Position)
# longevity_master$chromosomal_region = positions_longevity
# rsid_longevity = getBM(attributes = c('refsnp_id', 'chr_name',
#                                       'chrom_start'), 
#                         filters = 'chromosomal_region', 
#                         values = positions_longevity, 
#                         mart = snp_mart)
# rsid_longevity$chromosomal_region = paste0(rsid_longevity$chr_name, ":", rsid_longevity$chrom_start, ":", rsid_longevity$chrom_start)
# rsid_longevity_filtered = rsid_longevity[,c("refsnp_id", "chromosomal_region")]
# longevity_master$refsnp_id = rsid_longevity_filtered[match(longevity_master$chromosomal_region, rsid_longevity_filtered$chromosomal_region), "refsnp_id"]
# longevity = longevity_master[,c("refsnp_id", "NEA", "EA", "Effective_N", "Beta", "P.value")]

#low_density_lipoprotien: n/a
low_density_lipoprotien = low_density_lipoprotien_master[,c("rsID", "REF", "ALT", "N", "EFFECT_SIZE", "pvalue", "SE")]

#overweight: n
overweight_master$n = overweight_master$N_cases + overweight_master$N_controls
overweight = overweight_master[,c("MarkerName", "Allele2", "Allele1", "n", "beta", "P", "se")]

#percieved_age: n/a
percieved_age = percieved_age_master[,c("SNP", "ALLELE0", "ALLELE1", "N_EFFECTIVE", "LOG_OR", "P_BOLT_LMM_INF", "LOG_SE")]

#serum_urate_concentration: n
serum_urate_concentration_master$n = rep(110347, nrow(serum_urate_concentration_master))
serum_urate_concentration = serum_urate_concentration_master[,c("variant_id", "other_allele", "effect_allele", "n", "beta", "p_value", "standard_error")]

#telomere_length: n
telomere_length_master$n = rep(472174, nrow(telomere_length_master))
telomere_length = telomere_length_master[,c("variant_id", "other_allele", "effec_allele", "n", "beta", "p_value", "standard_error")]

#triglyceride: n/a
triglyceride = triglyceride_master[,c("rsID", "REF", "ALT", "N", "EFFECT_SIZE", "pvalue", "SE")]

#two_hour_glucose: n
two_hour_glucose_master$n = rep(15234, nrow(two_hour_glucose_master))
two_hour_glucose = two_hour_glucose_master[,c("snp", "other_allele", "effect_allele", "n", "effect", "pvalue", "stderr")]

#type_2_diabetes: n
type_2_diabetes_master$n = type_2_diabetes_master$N_CASES + type_2_diabetes_master$N_CONTROLS
type_2_diabetes = type_2_diabetes_master[,c("SNP", "OTHER_ALLELE", "RISK_ALLELE", "n", "OR", "P_VALUE")]
type_2_diabetes$standard_error = (type_2_diabetes_master$OR - type_2_diabetes_master$OR_95L)/1.96

#waist_hip_ratio: rsid
waist_hip_ratio_master$rsid = sub("\\:.*", "", waist_hip_ratio_master$SNP)
waist_hip_ratio = waist_hip_ratio_master[,c("rsid", "Other_Allele", "Tested_Allele", "N", "BETA", "P", "SE")]


#atherosclerosis: n/a
atherosclerosis = atherosclerosis_master[,c("variant_id", "other_allele", "effect_allele", "N", "beta", "p_value", "standard_error")]

#leucine: n/a
leucine = leucine_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#isoleucine: n/a
isoleucine = isoleucine_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#valine: n/a
valine = valine_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#free_cholesterol_in_very_large_hdl: n/a
free_cholesterol_in_very_large_hdl = free_cholesterol_in_very_large_hdl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#free_cholesterol_in_large_hdl: n/a
free_cholesterol_in_large_hdl = free_cholesterol_in_large_hdl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#phospholipid_in_very_large_hdl: n/a
phospholipid_in_very_large_hdl = phospholipid_in_very_large_hdl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#phospholipid_in_large_hdl: n/a
phospholipid_in_large_hdl = phospholipid_in_large_hdl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#triglyceride_in_very_large_vldl: n/a
triglyceride_in_very_large_vldl = triglyceride_in_very_large_vldl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]

#triglyceride_in_large_vldl: n/a
triglyceride_in_large_vldl = triglyceride_in_large_vldl_master[,c("ID", "NEA", "EA", "n_samples", "beta", "p.value", "se")]



print("2. calculate and filter all columns")

#==================
#3. rename all columns
#==================
new_column_names_ldsc_beta = c("SNP", "A2", "A1",  "N", "BETA", "P", "SE")
new_column_names_ldsc_or = c("SNP", "A2", "A1",  "N", "OR", "P", "SE")
new_column_names_ldsc_log_or = c("SNP", "A2", "A1",  "N", "LOG_ODDS", "P", "SE")

colnames(alzheimers) = new_column_names_ldsc_beta
colnames(body_mass_index) = new_column_names_ldsc_beta
colnames(chronic_inflammation) = new_column_names_ldsc_beta
colnames(class_1_obesity) = new_column_names_ldsc_beta
colnames(class_2_obesity) = new_column_names_ldsc_beta
colnames(class_3_obesity) = new_column_names_ldsc_beta
colnames(coronary_artery_disease) = new_column_names_ldsc_beta
colnames(extreme_body_mass_index) = new_column_names_ldsc_beta
colnames(extreme_waist_hip_ratio) = new_column_names_ldsc_beta
colnames(extrinsic_epigenetic_age_acceleration) = new_column_names_ldsc_beta
colnames(fasting_glucose) = new_column_names_ldsc_beta
colnames(fasting_insulin) = new_column_names_ldsc_beta
colnames(fasting_proinsulin) = new_column_names_ldsc_beta
colnames(frailty_index) = new_column_names_ldsc_beta
colnames(hba1c) = new_column_names_ldsc_beta
colnames(healthspan) = new_column_names_ldsc_beta
colnames(high_density_lipoprotien) = new_column_names_ldsc_beta
colnames(homab) = new_column_names_ldsc_beta
colnames(homair) = new_column_names_ldsc_beta
colnames(income) = new_column_names_ldsc_beta
colnames(intelligence) = new_column_names_ldsc_beta
colnames(intrinsic_epigenetic_age_acceleration) = new_column_names_ldsc_beta
colnames(lacunar_stroke) = new_column_names_ldsc_beta
colnames(low_density_lipoprotien) = new_column_names_ldsc_beta
colnames(overweight) = new_column_names_ldsc_beta
colnames(percieved_age) = new_column_names_ldsc_log_or #note: log odds ratio
colnames(serum_urate_concentration) = new_column_names_ldsc_beta
colnames(telomere_length) = new_column_names_ldsc_beta
colnames(triglyceride) = new_column_names_ldsc_beta
colnames(two_hour_glucose) = new_column_names_ldsc_beta
colnames(type_2_diabetes) = new_column_names_ldsc_or #note: odds ratio
colnames(waist_hip_ratio) = new_column_names_ldsc_beta
colnames(atherosclerosis) = new_column_names_ldsc_beta
colnames(leucine) = new_column_names_ldsc_beta
colnames(isoleucine) = new_column_names_ldsc_beta
colnames(valine) = new_column_names_ldsc_beta
colnames(free_cholesterol_in_very_large_hdl) = new_column_names_ldsc_beta
colnames(free_cholesterol_in_large_hdl) = new_column_names_ldsc_beta
colnames(phospholipid_in_very_large_hdl) = new_column_names_ldsc_beta
colnames(phospholipid_in_large_hdl) = new_column_names_ldsc_beta
colnames(triglyceride_in_very_large_vldl) = new_column_names_ldsc_beta
colnames(triglyceride_in_large_vldl) = new_column_names_ldsc_beta


print("3. rename all columns")

#==================
#4. filter rare and common snps
#==================

# #only need to filter one dataset for MAF, all datasets will be intersected afterwards
# fasting_glucose = fasting_glucose[which(fasting_glucose_master$maf >= 0.05),]
# 
# #find common snps and filter both datasets
# common_snps = Reduce(intersect, list(alzheimers$SNP, 
#                                      body_mass_index$SNP, 
#                                      chronic_inflammation$SNP, 
#                                      class_1_obesity$SNP, 
#                                      class_2_obesity$SNP, 
#                                      class_3_obesity$SNP, 
#                                      coronary_artery_disease$SNP, 
#                                      extreme_body_mass_index$SNP, 
#                                      extreme_waist_hip_ratio$SNP, 
#                                      extrinsic_epigenetic_age_acceleration$SNP, 
#                                      fasting_glucose$SNP, 
#                                      fasting_insulin$SNP, 
#                                      fasting_proinsulin$SNP, 
#                                      frailty_index$SNP, 
#                                      hba1c$SNP, 
#                                      healthspan$SNP, 
#                                      high_density_lipoprotien$SNP, 
#                                      homab$SNP, 
#                                      homair$SNP, 
#                                      income$SNP, 
#                                      intelligence$SNP, 
#                                      intrinsic_epigenetic_age_acceleration$SNP, 
#                                      lacunar_stroke$SNP,
#                                      low_density_lipoprotien$SNP, 
#                                      overweight$SNP, 
#                                      percieved_age$SNP, 
#                                      serum_urate_concentration$SNP, 
#                                      telomere_length$SNP, 
#                                      triglyceride$SNP, 
#                                      two_hour_glucose$SNP, 
#                                      type_2_diabetes$SNP, 
#                                      waist_hip_ratio$SNP,
#                                      atherosclerosis$SNP,
#                                      leucine$SNP, 
#                                      isoleucine$SNP,
#                                      valine$SNP, 
#                                      free_cholesterol_in_very_large_hdl$SNP, 
#                                      free_cholesterol_in_large_hdl$SNP, 
#                                      phospholipid_in_very_large_hdl$SNP, 
#                                      phospholipid_in_large_hdl$SNP, 
#                                      triglyceride_in_very_large_vldl$SNP, 
#                                      triglyceride_in_large_vldl$SNP))
# 
# 
# alzheimers = alzheimers[which(alzheimers$SNP %in% common_snps),]
# body_mass_index = body_mass_index[which(body_mass_index$SNP %in% common_snps),]
# chronic_inflammation = chronic_inflammation[which(chronic_inflammation$SNP %in% common_snps),]
# class_1_obesity = class_1_obesity[which(class_1_obesity$SNP %in% common_snps),]
# class_2_obesity = class_2_obesity[which(class_2_obesity$SNP %in% common_snps),]
# class_3_obesity = class_3_obesity[which(class_3_obesity$SNP %in% common_snps),]
# coronary_artery_disease = coronary_artery_disease[which(coronary_artery_disease$SNP %in% common_snps),]
# extreme_body_mass_index = extreme_body_mass_index[which(extreme_body_mass_index$SNP %in% common_snps),]
# extreme_waist_hip_ratio = extreme_waist_hip_ratio[which(extreme_waist_hip_ratio$SNP %in% common_snps),]
# extrinsic_epigenetic_age_acceleration = extrinsic_epigenetic_age_acceleration[which(extrinsic_epigenetic_age_acceleration$SNP %in% common_snps),]
# fasting_glucose = fasting_glucose[which(fasting_glucose$SNP %in% common_snps),]
# fasting_insulin = fasting_insulin[which(fasting_insulin$SNP %in% common_snps),]
# fasting_proinsulin = fasting_proinsulin[which(fasting_proinsulin$SNP %in% common_snps),]
# frailty_index = frailty_index[which(frailty_index$SNP %in% common_snps),]
# hba1c = hba1c[which(hba1c$SNP %in% common_snps),]
# healthspan = healthspan[which(healthspan$SNP %in% common_snps),]
# high_density_lipoprotien = high_density_lipoprotien[which(high_density_lipoprotien$SNP %in% common_snps),]
# homab = homab[which(homab$SNP %in% common_snps),]
# homair = homair[which(homair$SNP %in% common_snps),]
# income = income[which(income$SNP %in% common_snps),]
# intelligence = intelligence[which(intelligence$SNP %in% common_snps),]
# intrinsic_epigenetic_age_acceleration = intrinsic_epigenetic_age_acceleration[which(intrinsic_epigenetic_age_acceleration$SNP %in% common_snps),]
# lacunar_stroke = lacunar_stroke[which(lacunar_stroke$SNP %in% common_snps),]
# low_density_lipoprotien = low_density_lipoprotien[which(low_density_lipoprotien$SNP %in% common_snps),]
# overweight = overweight[which(overweight$SNP %in% common_snps),]
# percieved_age = percieved_age[which(percieved_age$SNP %in% common_snps),]
# serum_urate_concentration = serum_urate_concentration[which(serum_urate_concentration$SNP %in% common_snps),]
# telomere_length = telomere_length[which(telomere_length$SNP %in% common_snps),]
# triglyceride = triglyceride[which(triglyceride$SNP %in% common_snps),]
# two_hour_glucose = two_hour_glucose[which(two_hour_glucose$SNP %in% common_snps),]
# type_2_diabetes = type_2_diabetes[which(type_2_diabetes$SNP %in% common_snps),]
# waist_hip_ratio = waist_hip_ratio[which(waist_hip_ratio$SNP %in% common_snps),]
# atherosclerosis = atherosclerosis[which(atherosclerosis$SNP %in% common_snps),]
# leucine = leucine[which(leucine$SNP %in% common_snps),]
# isoleucine = isoleucine[which(isoleucine$SNP %in% common_snps),]
# valine = valine[which(valine$SNP %in% common_snps),]
# free_cholesterol_in_very_large_hdl = free_cholesterol_in_very_large_hdl[which(free_cholesterol_in_very_large_hdl$SNP %in% common_snps),]
# free_cholesterol_in_large_hdl = free_cholesterol_in_large_hdl[which(free_cholesterol_in_large_hdl$SNP %in% common_snps),]
# phospholipid_in_very_large_hdl = phospholipid_in_very_large_hdl[which(phospholipid_in_very_large_hdl$SNP %in% common_snps),]
# phospholipid_in_large_hdl = phospholipid_in_large_hdl[which(phospholipid_in_large_hdl$SNP %in% common_snps),]
# triglyceride_in_very_large_vldl = triglyceride_in_very_large_vldl[which(triglyceride_in_very_large_vldl$SNP %in% common_snps),]
# triglyceride_in_large_vldl = triglyceride_in_large_vldl[which(triglyceride_in_large_vldl$SNP %in% common_snps),]



alzheimers = alzheimers[!duplicated(alzheimers[,'SNP']),]
body_mass_index = body_mass_index[!duplicated(body_mass_index[,'SNP']),]
chronic_inflammation = chronic_inflammation[!duplicated(chronic_inflammation[,'SNP']),]
class_1_obesity = class_1_obesity[!duplicated(class_1_obesity[,'SNP']),]
class_2_obesity = class_2_obesity[!duplicated(class_2_obesity[,'SNP']),]
class_3_obesity = class_3_obesity[!duplicated(class_3_obesity[,'SNP']),]
coronary_artery_disease = coronary_artery_disease[!duplicated(coronary_artery_disease[,'SNP']),]
extreme_body_mass_index = extreme_body_mass_index[!duplicated(extreme_body_mass_index[,'SNP']),]
extreme_waist_hip_ratio = extreme_waist_hip_ratio[!duplicated(extreme_waist_hip_ratio[,'SNP']),]
extrinsic_epigenetic_age_acceleration = extrinsic_epigenetic_age_acceleration[!duplicated(extrinsic_epigenetic_age_acceleration[,'SNP']),]
fasting_glucose = fasting_glucose[!duplicated(fasting_glucose[,'SNP']),]
fasting_insulin = fasting_insulin[!duplicated(fasting_insulin[,'SNP']),]
fasting_proinsulin = fasting_proinsulin[!duplicated(fasting_proinsulin[,'SNP']),]
frailty_index = frailty_index[!duplicated(frailty_index[,'SNP']),]
hba1c = hba1c[!duplicated(hba1c[,'SNP']),]
healthspan = healthspan[!duplicated(healthspan[,'SNP']),]
high_density_lipoprotien = high_density_lipoprotien[!duplicated(high_density_lipoprotien[,'SNP']),]
homab = homab[!duplicated(homab[,'SNP']),]
homair = homair[!duplicated(homair[,'SNP']),]
income = income[!duplicated(income[,'SNP']),]
intelligence = intelligence[!duplicated(intelligence),]
intrinsic_epigenetic_age_acceleration = intrinsic_epigenetic_age_acceleration[!duplicated(intrinsic_epigenetic_age_acceleration[,'SNP']),]
lacunar_stroke = lacunar_stroke[!duplicated(lacunar_stroke[,'SNP']),]
low_density_lipoprotien = low_density_lipoprotien[!duplicated(low_density_lipoprotien[,'SNP']),]
overweight = overweight[!duplicated(overweight[,'SNP']),]
percieved_age = percieved_age[!duplicated(percieved_age[,'SNP']),]
serum_urate_concentration = serum_urate_concentration[!duplicated(serum_urate_concentration[,'SNP']),]
telomere_length = telomere_length[!duplicated(telomere_length[,'SNP']),]
triglyceride = triglyceride[!duplicated(triglyceride[,'SNP']),]
two_hour_glucose = two_hour_glucose[!duplicated(two_hour_glucose[,'SNP']),]
type_2_diabetes = type_2_diabetes[!duplicated(type_2_diabetes[,'SNP']),]
waist_hip_ratio = waist_hip_ratio[!duplicated(waist_hip_ratio[,'SNP']),]
atherosclerosis = atherosclerosis[!duplicated(atherosclerosis[,'SNP']),]
leucine = leucine[!duplicated(leucine[,'SNP']),]
isoleucine = isoleucine[!duplicated(isoleucine[,'SNP']),]
valine = valine[!duplicated(valine[,'SNP']),]
free_cholesterol_in_very_large_hdl = free_cholesterol_in_very_large_hdl[!duplicated(free_cholesterol_in_very_large_hdl[,'SNP']),]
free_cholesterol_in_large_hdl = free_cholesterol_in_large_hdl[!duplicated(free_cholesterol_in_large_hdl[,'SNP']),]
phospholipid_in_very_large_hdl = phospholipid_in_very_large_hdl[!duplicated(phospholipid_in_very_large_hdl[,'SNP']),]
phospholipid_in_large_hdl = phospholipid_in_large_hdl[!duplicated(phospholipid_in_large_hdl[,'SNP']),]
triglyceride_in_very_large_vldl = triglyceride_in_very_large_vldl[!duplicated(triglyceride_in_very_large_vldl[,'SNP']),]
triglyceride_in_large_vldl = triglyceride_in_large_vldl[!duplicated(triglyceride_in_large_vldl[,'SNP']),]


print("4. filter rare and common snps")

#==================
#5. save all files
#==================

# homab$BETA = -1 * homab$BETA

print(max(alzheimers$N)) 
print(max(body_mass_index$N)) 
print(max(chronic_inflammation$N)) 
print(max(class_1_obesity$N)) 
print(max(class_2_obesity$N)) 
print(max(class_3_obesity$N)) 
print(max(coronary_artery_disease$N)) 
print(max(extreme_body_mass_index$N)) 
print(max(extreme_waist_hip_ratio$N)) 
print(max(extrinsic_epigenetic_age_acceleration$N)) 
print(max(fasting_glucose$N)) 
print(max(fasting_insulin$N)) 
print(max(fasting_proinsulin$N)) 
print(max(frailty_index$N)) 
print(max(hba1c$N)) 
print(max(healthspan$N)) 
print(max(high_density_lipoprotien$N)) 
print(max(homab$N)) 
print(max(homair$N)) 
print(max(income$N)) 
print(max(intelligence$N)) 
print(max(intrinsic_epigenetic_age_acceleration$N)) 
print(max(lacunar_stroke$N))
print(max(low_density_lipoprotien$N)) 
print(max(overweight$N)) 
print(max(percieved_age$N)) 
print(max(serum_urate_concentration$N)) 
print(max(telomere_length$N)) 
print(max(triglyceride$N)) 
print(max(two_hour_glucose$N)) 
print(max(type_2_diabetes$N)) 
print(max(waist_hip_ratio$N))
print(max(atherosclerosis$N))
print(max(leucine$N)) 
print(max(isoleucine$N))
print(max(valine$N)) 
print(max(free_cholesterol_in_very_large_hdl$N)) 
print(max(free_cholesterol_in_large_hdl$N)) 
print(max(phospholipid_in_very_large_hdl$N)) 
print(max(phospholipid_in_large_hdl$N)) 
print(max(triglyceride_in_very_large_vldl$N)) 
print(max(triglyceride_in_large_vldl$N))



write.table(alzheimers,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/alzheimers.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(body_mass_index,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/body_mass_index.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(chronic_inflammation,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/chronic_inflammation.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(class_1_obesity,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_1_obesity.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(class_2_obesity,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_2_obesity.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(class_3_obesity,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_3_obesity.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(coronary_artery_disease,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/coronary_artery_disease.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(extreme_body_mass_index,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_body_mass_index.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(extreme_waist_hip_ratio,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_waist_hip_ratio.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(extrinsic_epigenetic_age_acceleration,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extrinsic_epigenetic_age_acceleration.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(fasting_glucose,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_glucose.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(fasting_insulin,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_insulin.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(fasting_proinsulin,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_proinsulin.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(frailty_index,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/frailty_index.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(hba1c,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/hba1c.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(healthspan,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/healthspan.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(high_density_lipoprotien,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/high_density_lipoprotien.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(homab,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homab.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(homair,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homair.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(income,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/income.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(intelligence,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intelligence.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(intrinsic_epigenetic_age_acceleration,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intrinsic_epigenetic_age_acceleration.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(lacunar_stroke,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/lacunar_stroke.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(low_density_lipoprotien,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/low_density_lipoprotien.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(overweight,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/overweight.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(percieved_age,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/percieved_age.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(serum_urate_concentration,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/serum_urate_concentration.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(telomere_length,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/telomere_length.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(triglyceride,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(two_hour_glucose,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/two_hour_glucose.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(type_2_diabetes,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/type_2_diabetes.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(waist_hip_ratio,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/waist_hip_ratio.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(atherosclerosis,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/atherosclerosis.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(leucine,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/leucine.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(isoleucine,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/isoleucine.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(valine,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/valine.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(free_cholesterol_in_very_large_hdl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_very_large_hdl.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(free_cholesterol_in_large_hdl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_large_hdl.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(phospholipid_in_very_large_hdl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_very_large_hdl.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(phospholipid_in_large_hdl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_large_hdl.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(triglyceride_in_very_large_vldl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_very_large_vldl.txt',sep=" ",row.names=FALSE, quote=FALSE)
write.table(triglyceride_in_large_vldl,'~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_large_vldl.txt',sep=" ",row.names=FALSE, quote=FALSE)


print("5. save all files")
