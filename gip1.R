

alzheimers = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/alzheimers.txt', header = T)
body_mass_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/body_mass_index.txt', header = T)
chronic_inflammation = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/chronic_inflammation.txt', header = T)
class_1_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_1_obesity.txt', header = T)
class_2_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_2_obesity.txt', header = T)
class_3_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_3_obesity.txt', header = T)
coronary_artery_disease = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/coronary_artery_disease.txt', header = T)
extreme_body_mass_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_body_mass_index.txt', header = T)
extreme_waist_hip_ratio = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_waist_hip_ratio.txt', header = T)
extrinsic_epigenetic_age_acceleration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extrinsic_epigenetic_age_acceleration.txt', header = T)
fasting_glucose = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_glucose.txt', header = T)
fasting_insulin = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_insulin.txt', header = T)
fasting_proinsulin = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_proinsulin.txt', header = T)
frailty_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/frailty_index.txt', header = T)
hba1c = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/hba1c.txt', header = T)
healthspan = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/healthspan.txt', header = T)
high_density_lipoprotien = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/high_density_lipoprotien.txt', header = T)
homab = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homab.txt', header = T)
homair = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homair.txt', header = T)
income = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/income.txt', header = T)
intelligence = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intelligence.txt', header = T)
intrinsic_epigenetic_age_acceleration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intrinsic_epigenetic_age_acceleration.txt', header = T)
lacunar_stroke = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/lacunar_stroke.txt', header = T)
low_density_lipoprotien = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/low_density_lipoprotien.txt', header = T)
overweight = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/overweight.txt', header = T)
percieved_age = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/percieved_age.txt', header = T)
serum_urate_concentration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/serum_urate_concentration.txt', header = T)
telomere_length = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/telomere_length.txt', header = T)
triglyceride = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride.txt', header = T)
two_hour_glucose = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/two_hour_glucose.txt', header = T)
type_2_diabetes = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/type_2_diabetes.txt', header = T)
waist_hip_ratio = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/waist_hip_ratio.txt', header = T)
atherosclerosis = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/atherosclerosis.txt', header = T)
leucine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/leucine.txt', header = T)
isoleucine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/isoleucine.txt', header = T)
valine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/valine.txt', header = T)
free_cholesterol_in_very_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_very_large_hdl.txt', header = T)
free_cholesterol_in_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_large_hdl.txt', header = T)
phospholipid_in_very_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_very_large_hdl.txt', header = T)
phospholipid_in_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_large_hdl.txt', header = T)
triglyceride_in_very_large_vldl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_very_large_vldl.txt', header = T)
triglyceride_in_large_vldl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_large_vldl.txt', header = T)




# alzheimers = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/alzheimers.txt', header = T, nrow=20)
# body_mass_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/body_mass_index.txt', header = T, nrow=20)
# chronic_inflammation = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/chronic_inflammation.txt', header = T, nrow=20)
# class_1_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_1_obesity.txt', header = T, nrow=20)
# class_2_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_2_obesity.txt', header = T, nrow=20)
# class_3_obesity = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/class_3_obesity.txt', header = T, nrow=20)
# coronary_artery_disease = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/coronary_artery_disease.txt', header = T, nrow=20)
# extreme_body_mass_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_body_mass_index.txt', header = T, nrow=20)
# extreme_waist_hip_ratio = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extreme_waist_hip_ratio.txt', header = T, nrow=20)
# extrinsic_epigenetic_age_acceleration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/extrinsic_epigenetic_age_acceleration.txt', header = T, nrow=20)
# fasting_glucose = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_glucose.txt', header = T, nrow=20)
# fasting_insulin = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_insulin.txt', header = T, nrow=20)
# fasting_proinsulin = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/fasting_proinsulin.txt', header = T, nrow=20)
# frailty_index = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/frailty_index.txt', header = T, nrow=20)
# hba1c = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/hba1c.txt', header = T, nrow=20)
# healthspan = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/healthspan.txt', header = T, nrow=20)
# high_density_lipoprotien = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/high_density_lipoprotien.txt', header = T, nrow=20)
# homab = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homab.txt', header = T, nrow=20)
# homair = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/homair.txt', header = T, nrow=20)
# income = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/income.txt', header = T, nrow=20)
# intelligence = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intelligence.txt', header = T, nrow=20)
# intrinsic_epigenetic_age_acceleration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/intrinsic_epigenetic_age_acceleration.txt', header = T, nrow=20)
# lacunar_stroke = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/lacunar_stroke.txt', header = T, nrow=20)
# low_density_lipoprotien = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/low_density_lipoprotien.txt', header = T, nrow=20)
# overweight = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/overweight.txt', header = T, nrow=20)
# percieved_age = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/percieved_age.txt', header = T, nrow=20)
# serum_urate_concentration = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/serum_urate_concentration.txt', header = T, nrow=20)
# telomere_length = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/telomere_length.txt', header = T, nrow=20)
# triglyceride = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride.txt', header = T, nrow=20)
# two_hour_glucose = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/two_hour_glucose.txt', header = T, nrow=20)
# type_2_diabetes = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/type_2_diabetes.txt', header = T, nrow=20)
# waist_hip_ratio = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/waist_hip_ratio.txt', header = T, nrow=20)
# atherosclerosis = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/atherosclerosis.txt', header = T, nrow=20)
# leucine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/leucine.txt', header = T, nrow=20)
# isoleucine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/isoleucine.txt', header = T, nrow=20)
# valine = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/valine.txt', header = T, nrow=20)
# free_cholesterol_in_very_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_very_large_hdl.txt', header = T, nrow=20)
# free_cholesterol_in_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/free_cholesterol_in_large_hdl.txt', header = T, nrow=20)
# phospholipid_in_very_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_very_large_hdl.txt', header = T, nrow=20)
# phospholipid_in_large_hdl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/phospholipid_in_large_hdl.txt', header = T, nrow=20)
# triglyceride_in_very_large_vldl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_very_large_vldl.txt', header = T, nrow=20)
# triglyceride_in_large_vldl = read.table('~/scratch/principal_component_GWAS/preprocessed_with_standard_error/triglyceride_in_large_vldl.txt', header = T, nrow=20)



#find common snps and filter both datasets
common_snps = Reduce(intersect, list(#alzheimers$SNP,
                                     body_mass_index$SNP,
                                     chronic_inflammation$SNP,
                                     #class_1_obesity$SNP,
                                     #class_2_obesity$SNP,
                                     #class_3_obesity$SNP,
                                     coronary_artery_disease$SNP,
                                     #extreme_body_mass_index$SNP,
                                     #extreme_waist_hip_ratio$SNP,
                                     #extrinsic_epigenetic_age_acceleration$SNP,
                                     fasting_glucose$SNP,
                                     fasting_insulin$SNP,
                                     fasting_proinsulin$SNP,
                                     #frailty_index$SNP,
                                     hba1c$SNP,
                                     healthspan$SNP,
                                     high_density_lipoprotien$SNP,
                                     homab$SNP,
                                     homair$SNP,
                                     #income$SNP,
                                     #intelligence$SNP,
                                     #intrinsic_epigenetic_age_acceleration$SNP,
                                     lacunar_stroke$SNP,
                                     low_density_lipoprotien$SNP,
                                     #overweight$SNP,
                                     #percieved_age$SNP,
                                     serum_urate_concentration$SNP,
                                     #telomere_length$SNP,
                                     triglyceride$SNP,
                                     two_hour_glucose$SNP,
                                     type_2_diabetes$SNP,
                                     waist_hip_ratio$SNP,
                                     atherosclerosis$SNP,
                                     leucine$SNP,
                                     isoleucine$SNP,
                                     valine$SNP#,
                                     #free_cholesterol_in_very_large_hdl$SNP,
                                     #free_cholesterol_in_large_hdl$SNP,
                                     #phospholipid_in_very_large_hdl$SNP,
                                     #phospholipid_in_large_hdl$SNP,
                                     #triglyceride_in_very_large_vldl$SNP,
                                     #triglyceride_in_large_vldl$SNP
                                     ))




# alzheimers = alzheimers[which(alzheimers$SNP %in% common_snps),]
body_mass_index = body_mass_index[which(body_mass_index$SNP %in% common_snps),]
chronic_inflammation = chronic_inflammation[which(chronic_inflammation$SNP %in% common_snps),]
# class_1_obesity = class_1_obesity[which(class_1_obesity$SNP %in% common_snps),]
# class_2_obesity = class_2_obesity[which(class_2_obesity$SNP %in% common_snps),]
# class_3_obesity = class_3_obesity[which(class_3_obesity$SNP %in% common_snps),]
coronary_artery_disease = coronary_artery_disease[which(coronary_artery_disease$SNP %in% common_snps),]
# extreme_body_mass_index = extreme_body_mass_index[which(extreme_body_mass_index$SNP %in% common_snps),]
# extreme_waist_hip_ratio = extreme_waist_hip_ratio[which(extreme_waist_hip_ratio$SNP %in% common_snps),]
# extrinsic_epigenetic_age_acceleration = extrinsic_epigenetic_age_acceleration[which(extrinsic_epigenetic_age_acceleration$SNP %in% common_snps),]
fasting_glucose = fasting_glucose[which(fasting_glucose$SNP %in% common_snps),]
fasting_insulin = fasting_insulin[which(fasting_insulin$SNP %in% common_snps),]
fasting_proinsulin = fasting_proinsulin[which(fasting_proinsulin$SNP %in% common_snps),]
# frailty_index = frailty_index[which(frailty_index$SNP %in% common_snps),]
hba1c = hba1c[which(hba1c$SNP %in% common_snps),]
healthspan = healthspan[which(healthspan$SNP %in% common_snps),]
high_density_lipoprotien = high_density_lipoprotien[which(high_density_lipoprotien$SNP %in% common_snps),]
homab = homab[which(homab$SNP %in% common_snps),]
homair = homair[which(homair$SNP %in% common_snps),]
# income = income[which(income$SNP %in% common_snps),]
# intelligence = intelligence[which(intelligence$SNP %in% common_snps),]
# intrinsic_epigenetic_age_acceleration = intrinsic_epigenetic_age_acceleration[which(intrinsic_epigenetic_age_acceleration$SNP %in% common_snps),]
lacunar_stroke = lacunar_stroke[which(lacunar_stroke$SNP %in% common_snps),]
low_density_lipoprotien = low_density_lipoprotien[which(low_density_lipoprotien$SNP %in% common_snps),]
# overweight = overweight[which(overweight$SNP %in% common_snps),]
# percieved_age = percieved_age[which(percieved_age$SNP %in% common_snps),]
serum_urate_concentration = serum_urate_concentration[which(serum_urate_concentration$SNP %in% common_snps),]
# telomere_length = telomere_length[which(telomere_length$SNP %in% common_snps),]
triglyceride = triglyceride[which(triglyceride$SNP %in% common_snps),]
two_hour_glucose = two_hour_glucose[which(two_hour_glucose$SNP %in% common_snps),]
type_2_diabetes = type_2_diabetes[which(type_2_diabetes$SNP %in% common_snps),]
waist_hip_ratio = waist_hip_ratio[which(waist_hip_ratio$SNP %in% common_snps),]
atherosclerosis = atherosclerosis[which(atherosclerosis$SNP %in% common_snps),]
leucine = leucine[which(leucine$SNP %in% common_snps),]
isoleucine = isoleucine[which(isoleucine$SNP %in% common_snps),]
valine = valine[which(valine$SNP %in% common_snps),]
# free_cholesterol_in_very_large_hdl = free_cholesterol_in_very_large_hdl[which(free_cholesterol_in_very_large_hdl$SNP %in% common_snps),]
# free_cholesterol_in_large_hdl = free_cholesterol_in_large_hdl[which(free_cholesterol_in_large_hdl$SNP %in% common_snps),]
# phospholipid_in_very_large_hdl = phospholipid_in_very_large_hdl[which(phospholipid_in_very_large_hdl$SNP %in% common_snps),]
# phospholipid_in_large_hdl = phospholipid_in_large_hdl[which(phospholipid_in_large_hdl$SNP %in% common_snps),]
# triglyceride_in_very_large_vldl = triglyceride_in_very_large_vldl[which(triglyceride_in_very_large_vldl$SNP %in% common_snps),]
# triglyceride_in_large_vldl = triglyceride_in_large_vldl[which(triglyceride_in_large_vldl$SNP %in% common_snps),]


print(1)


#1. Estimated Ω and ΩSE using LDSC software
ldsc_results = read.table("~/scratch/principal_component_GWAS/all_traits_raw/results.txt", header = T)
library(stringr)
ldsc_results$trait_1 = unlist(lapply(ldsc_results$p1, function(x) str_match(x, "munge/munge_out/\\s*(.*?)\\s*.sumstats.gz")[,2]))
ldsc_results$trait_2 = unlist(lapply(ldsc_results$p2, function(x) str_match(x, "munge/munge_out/\\s*(.*?)\\s*.sumstats.gz")[,2]))





traits_to_keep = c("homair",
                   "homab",
                   "hba1c",
                   "fasting_glucose",
                   "fasting_insulin",
                   "fasting_proinsulin",
                   "triglyceride",
                   "serum_urate_concentration",
                   "waist_hip_ratio",
                   "body_mass_index",
                   "chronic_inflammation",
                   "valine",
                   "type_2_diabetes",
                   "leucine",
                   "isoleucine",
                   "healthspan",
                   "atherosclerosis",
                   "coronary_artery_disease",
                   "low_density_lipoprotien",
                   "high_density_lipoprotien",
                   "lacunar_stroke",
                   "two_hour_glucose"
                   )

ldsc_results = ldsc_results[which((ldsc_results$trait_1 %in% traits_to_keep)),]
ldsc_results = ldsc_results[which((ldsc_results$trait_2 %in% traits_to_keep)),]



library(maditr)
temp = ldsc_results[,c("trait_1", "trait_2", "rg")]
temp_opposite = temp
colnames(temp_opposite) = c("trait_2", "trait_1", "rg")
temp_opposite = temp_opposite[,c("trait_1", "trait_2", "rg")]



temp_full = rbind(temp, temp_opposite)
casted = dcast(temp_full, trait_1 ~ trait_2, value.var = "rg")

temp_correlation = as.matrix(casted)
rownames(temp_correlation) = temp_correlation[,1]
temp_correlation = temp_correlation[,-1]
temp_correlation = `class<-`(temp_correlation, 'numeric')


diag(temp_correlation) = 1


temp_correlation[which(is.na(as.numeric(as.character(temp_correlation))))] = 0  

#reorder correlation amtrix to our beta values order
temp_correlation = temp_correlation[traits_to_keep, traits_to_keep]
# 
# #===== calculate standard error
# 
# 
# library(maditr)
# temp = ldsc_results[,c("trait_1", "trait_2", "se")]
# temp_opposite = temp
# colnames(temp_opposite) = c("trait_2", "trait_1", "se")
# temp_opposite = temp_opposite[,c("trait_1", "trait_2", "se")]
# 
# 
# 
# temp_full = rbind(temp, temp_opposite)
# casted = dcast(temp_full, trait_1 ~ trait_2, value.var = "se")
# 
# temp_standard_error = as.matrix(casted)
# rownames(temp_standard_error) = temp_standard_error[,1]
# temp_standard_error = temp_standard_error[,-1]
# temp_standard_error = `class<-`(temp_standard_error, 'numeric')
# 
# 
# diag(temp_standard_error) = 1
# 
# temp_standard_error[which(is.na(as.numeric(as.character(temp_standard_error))))] = 0  
# 

print(2)

#2. Estimated varYi and Pearson correlation matrix for four pain phenotypes.


print(3)

# #3. Standardized GWAS summary statistics for four pain phenotypes (𝛽𝑠 = 𝛽 ⁄𝑆𝐷 and 𝑆𝐸𝑠 = 𝑆𝐸 ⁄𝑆𝐷 ).

#create beta and se for odds ratio and CI traits
#note: odds ratio (type 2 diabetes)
type_2_diabetes$BETA = log(type_2_diabetes$OR)
#note: log odds ratio (percieved age)
percieved_age$OR = exp(percieved_age$LOG_ODDS)
percieved_age$BETA = log(percieved_age$OR)

  
#=====
#scale start
#====

# 
# # alzheimers$BETA_STANDARDIZED = scale(alzheimers$BETA) 
# body_mass_index$BETA_STANDARDIZED = scale(body_mass_index$BETA) 
# chronic_inflammation$BETA_STANDARDIZED = scale(chronic_inflammation$BETA) 
# # class_1_obesity$BETA_STANDARDIZED = scale(class_1_obesity$BETA) 
# # class_2_obesity$BETA_STANDARDIZED = scale(class_2_obesity$BETA) 
# # class_3_obesity$BETA_STANDARDIZED = scale(class_3_obesity$BETA) 
# coronary_artery_disease$BETA_STANDARDIZED = scale(coronary_artery_disease$BETA) 
# # extreme_body_mass_index$BETA_STANDARDIZED = scale(extreme_body_mass_index$BETA) 
# # extreme_waist_hip_ratio$BETA_STANDARDIZED = scale(extreme_waist_hip_ratio$BETA) 
# # extrinsic_epigenetic_age_acceleration$BETA_STANDARDIZED = scale(extrinsic_epigenetic_age_acceleration$BETA) 
# fasting_glucose$BETA_STANDARDIZED = scale(fasting_glucose$BETA) 
# fasting_insulin$BETA_STANDARDIZED = scale(fasting_insulin$BETA) 
# fasting_proinsulin$BETA_STANDARDIZED = scale(fasting_proinsulin$BETA) 
# # frailty_index$BETA_STANDARDIZED = scale(frailty_index$BETA) 
# hba1c$BETA_STANDARDIZED = scale(hba1c$BETA) 
# healthspan$BETA_STANDARDIZED = scale(healthspan$BETA) 
# high_density_lipoprotien$BETA_STANDARDIZED = scale(high_density_lipoprotien$BETA) 
# homab$BETA_STANDARDIZED = scale(homab$BETA) 
# homair$BETA_STANDARDIZED = scale(homair$BETA) 
# # income$BETA_STANDARDIZED = scale(income$BETA) 
# # intelligence$BETA_STANDARDIZED = scale(intelligence$BETA) 
# # intrinsic_epigenetic_age_acceleration$BETA_STANDARDIZED = scale(intrinsic_epigenetic_age_acceleration$BETA) 
# lacunar_stroke$BETA_STANDARDIZED = scale(lacunar_stroke$BETA)
# low_density_lipoprotien$BETA_STANDARDIZED = scale(low_density_lipoprotien$BETA) 
# # overweight$BETA_STANDARDIZED = scale(overweight$BETA) 
# # percieved_age$BETA_STANDARDIZED = scale(percieved_age$BETA) 
# serum_urate_concentration$BETA_STANDARDIZED = scale(serum_urate_concentration$BETA) 
# # telomere_length$BETA_STANDARDIZED = scale(telomere_length$BETA) 
# triglyceride$BETA_STANDARDIZED = scale(triglyceride$BETA) 
# two_hour_glucose$BETA_STANDARDIZED = scale(two_hour_glucose$BETA) 
# type_2_diabetes$BETA_STANDARDIZED = scale(type_2_diabetes$BETA) 
# waist_hip_ratio$BETA_STANDARDIZED = scale(waist_hip_ratio$BETA)
# atherosclerosis$BETA_STANDARDIZED = scale(atherosclerosis$BETA)
# leucine$BETA_STANDARDIZED = scale(leucine$BETA) 
# isoleucine$BETA_STANDARDIZED = scale(isoleucine$BETA)
# valine$BETA_STANDARDIZED = scale(valine$BETA) 
# # free_cholesterol_in_very_large_hdl$BETA_STANDARDIZED = scale(free_cholesterol_in_very_large_hdl$BETA) 
# # free_cholesterol_in_large_hdl$BETA_STANDARDIZED = scale(free_cholesterol_in_large_hdl$BETA) 
# # phospholipid_in_very_large_hdl$BETA_STANDARDIZED = scale(phospholipid_in_very_large_hdl$BETA) 
# # phospholipid_in_large_hdl$BETA_STANDARDIZED = scale(phospholipid_in_large_hdl$BETA) 
# # triglyceride_in_very_large_vldl$BETA_STANDARDIZED = scale(triglyceride_in_very_large_vldl$BETA) 
# # triglyceride_in_large_vldl$BETA_STANDARDIZED = scale(triglyceride_in_large_vldl$BETA)
# 
# 
# 
# # alzheimers$SE_STANDARDIZED = scale(alzheimers$SE) 
# body_mass_index$SE_STANDARDIZED = scale(body_mass_index$SE) 
# chronic_inflammation$SE_STANDARDIZED = scale(chronic_inflammation$SE) 
# # class_1_obesity$SE_STANDARDIZED = scale(class_1_obesity$SE) 
# # class_2_obesity$SE_STANDARDIZED = scale(class_2_obesity$SE) 
# # class_3_obesity$SE_STANDARDIZED = scale(class_3_obesity$SE) 
# coronary_artery_disease$SE_STANDARDIZED = scale(coronary_artery_disease$SE) 
# # extreme_body_mass_index$SE_STANDARDIZED = scale(extreme_body_mass_index$SE) 
# # extreme_waist_hip_ratio$SE_STANDARDIZED = scale(extreme_waist_hip_ratio$SE) 
# # extrinsic_epigenetic_age_acceleration$SE_STANDARDIZED = scale(extrinsic_epigenetic_age_acceleration$SE) 
# fasting_glucose$SE_STANDARDIZED = scale(fasting_glucose$SE) 
# fasting_insulin$SE_STANDARDIZED = scale(fasting_insulin$SE) 
# fasting_proinsulin$SE_STANDARDIZED = scale(fasting_proinsulin$SE) 
# # frailty_index$SE_STANDARDIZED = scale(frailty_index$SE) 
# hba1c$SE_STANDARDIZED = scale(hba1c$SE) 
# healthspan$SE_STANDARDIZED = scale(healthspan$SE) 
# high_density_lipoprotien$SE_STANDARDIZED = scale(high_density_lipoprotien$SE) 
# homab$SE_STANDARDIZED = scale(homab$SE) 
# homair$SE_STANDARDIZED = scale(homair$SE) 
# # income$SE_STANDARDIZED = scale(income$SE) 
# # intelligence$SE_STANDARDIZED = scale(intelligence$SE) 
# # intrinsic_epigenetic_age_acceleration$SE_STANDARDIZED = scale(intrinsic_epigenetic_age_acceleration$SE) 
# lacunar_stroke$SE_STANDARDIZED = scale(lacunar_stroke$SE)
# low_density_lipoprotien$SE_STANDARDIZED = scale(low_density_lipoprotien$SE) 
# # overweight$SE_STANDARDIZED = scale(overweight$SE) 
# # percieved_age$SE_STANDARDIZED = scale(percieved_age$SE) 
# serum_urate_concentration$SE_STANDARDIZED = scale(serum_urate_concentration$SE) 
# # telomere_length$SE_STANDARDIZED = scale(telomere_length$SE) 
# triglyceride$SE_STANDARDIZED = scale(triglyceride$SE) 
# two_hour_glucose$SE_STANDARDIZED = scale(two_hour_glucose$SE) 
# type_2_diabetes$SE_STANDARDIZED = scale(type_2_diabetes$SE) 
# waist_hip_ratio$SE_STANDARDIZED = scale(waist_hip_ratio$SE)
# atherosclerosis$SE_STANDARDIZED = scale(atherosclerosis$SE)
# leucine$SE_STANDARDIZED = scale(leucine$SE) 
# isoleucine$SE_STANDARDIZED = scale(isoleucine$SE)
# valine$SE_STANDARDIZED = scale(valine$SE) 
# # free_cholesterol_in_very_large_hdl$SE_STANDARDIZED = scale(free_cholesterol_in_very_large_hdl$SE) 
# # free_cholesterol_in_large_hdl$SE_STANDARDIZED = scale(free_cholesterol_in_large_hdl$SE) 
# # phospholipid_in_very_large_hdl$SE_STANDARDIZED = scale(phospholipid_in_very_large_hdl$SE) 
# # phospholipid_in_large_hdl$SE_STANDARDIZED = scale(phospholipid_in_large_hdl$SE) 
# # triglyceride_in_very_large_vldl$SE_STANDARDIZED = scale(triglyceride_in_very_large_vldl$SE) 
# # triglyceride_in_large_vldl$SE_STANDARDIZED = scale(triglyceride_in_large_vldl$SE)

#=====
#scale end
#====







#4. Checked whether all eigenvalues were positive for Ω.
# print(eigenvalues)

#5. Estimated eigenvalues (L) and the matrix of eigenvectors (A) of Ω.

eigenvectors = prcomp(temp_correlation, scale = T)$rotation
eigenvalues = prcomp(temp_correlation, scale = T)$center


print(4)

#6. If the coefficient of a given eigenvector for back pain was negative (𝑎𝑖,𝑏𝑎𝑐𝑘 𝑝𝑎𝑖𝑛< 0), we changed the signs for all coefficients in its eigenvector
# eigenvectors = abs(eigenvectors)
eigenvectors[which(eigenvalues < 0),] = -1 * eigenvectors[which(eigenvalues < 0),]
eigenvalues[which(eigenvalues < 0)] = -1 * eigenvalues[which(eigenvalues < 0)]




#7. Estimated variance for GIP as 𝑣𝑎𝑟(𝐺𝑂𝑃 ) = ∑[(𝑎 ⨂𝑎 ) ∘ ∑], where ⨂ is an outer product.
matrix_temp = as.matrix((eigenvectors[,1] %o% eigenvectors[,1])*temp_correlation)
matrix_temp[is.na(matrix_temp)] = 0

gip1_variance = sum(matrix_temp)

#8. Scaled coefficients for GIPs as 𝑎𝑠 = 𝑎 ⁄𝑆𝐷(𝐺𝑂𝑃 ).

#change_me try without this step to see if log p values are smaller
#eigenvectors = eigenvectors/sqrt(gip1_variance)





print(5)

# 
# 
# beta_matrix = matrix(c(alzheimers$BETA_STANDARDIZED,
#                         body_mass_index$BETA_STANDARDIZED,
#                         chronic_inflammation$BETA_STANDARDIZED,
#                         class_1_obesity$BETA_STANDARDIZED,
#                         class_2_obesity$BETA_STANDARDIZED,
#                         class_3_obesity$BETA_STANDARDIZED,
#                         coronary_artery_disease$BETA_STANDARDIZED,
#                         extreme_body_mass_index$BETA_STANDARDIZED,
#                         extreme_waist_hip_ratio$BETA_STANDARDIZED,
#                         extrinsic_epigenetic_age_acceleration$BETA_STANDARDIZED,
#                         fasting_glucose$BETA_STANDARDIZED,
#                         fasting_insulin$BETA_STANDARDIZED,
#                         fasting_proinsulin$BETA_STANDARDIZED,
#                         frailty_index$BETA_STANDARDIZED,
#                         hba1c$BETA_STANDARDIZED,
#                         healthspan$BETA_STANDARDIZED,
#                         high_density_lipoprotien$BETA_STANDARDIZED,
#                         homab$BETA_STANDARDIZED,
#                         homair$BETA_STANDARDIZED,
#                         income$BETA_STANDARDIZED,
#                         intelligence$BETA_STANDARDIZED,
#                         intrinsic_epigenetic_age_acceleration$BETA_STANDARDIZED,
#                         lacunar_stroke$BETA_STANDARDIZED,
#                         low_density_lipoprotien$BETA_STANDARDIZED,
#                         overweight$BETA_STANDARDIZED,
#                         percieved_age$BETA_STANDARDIZED,
#                         serum_urate_concentration$BETA_STANDARDIZED,
#                         telomere_length$BETA_STANDARDIZED,
#                         triglyceride$BETA_STANDARDIZED,
#                         two_hour_glucose$BETA_STANDARDIZED,
#                         type_2_diabetes$BETA_STANDARDIZED,
#                         waist_hip_ratio$BETA_STANDARDIZED,
#                         atherosclerosis$BETA_STANDARDIZED,
#                         leucine$BETA_STANDARDIZED,
#                         isoleucine$BETA_STANDARDIZED,
#                         valine$BETA_STANDARDIZED,
#                         free_cholesterol_in_very_large_hdl$BETA_STANDARDIZED,
#                         free_cholesterol_in_large_hdl$BETA_STANDARDIZED,
#                         phospholipid_in_very_large_hdl$BETA_STANDARDIZED,
#                         phospholipid_in_large_hdl$BETA_STANDARDIZED,
#                         triglyceride_in_very_large_vldl$BETA_STANDARDIZED,
#                         triglyceride_in_large_vldl$BETA_STANDARDIZED), ncol = 42, byrow = F)
# 
# se_matrix = matrix(c(alzheimers$SE_STANDARDIZED,
#                      body_mass_index$SE_STANDARDIZED,
#                      chronic_inflammation$SE_STANDARDIZED,
#                      class_1_obesity$SE_STANDARDIZED,
#                      class_2_obesity$SE_STANDARDIZED,
#                      class_3_obesity$SE_STANDARDIZED,
#                      coronary_artery_disease$SE_STANDARDIZED,
#                      extreme_body_mass_index$SE_STANDARDIZED,
#                      extreme_waist_hip_ratio$SE_STANDARDIZED,
#                      extrinsic_epigenetic_age_acceleration$SE_STANDARDIZED,
#                      fasting_glucose$SE_STANDARDIZED,
#                      fasting_insulin$SE_STANDARDIZED,
#                      fasting_proinsulin$SE_STANDARDIZED,
#                      frailty_index$SE_STANDARDIZED,
#                      hba1c$SE_STANDARDIZED,
#                      healthspan$SE_STANDARDIZED,
#                      high_density_lipoprotien$SE_STANDARDIZED,
#                      homab$SE_STANDARDIZED,
#                      homair$SE_STANDARDIZED,
#                      income$SE_STANDARDIZED,
#                      intelligence$SE_STANDARDIZED,
#                      intrinsic_epigenetic_age_acceleration$SE_STANDARDIZED,
#                      lacunar_stroke$SE_STANDARDIZED,
#                      low_density_lipoprotien$SE_STANDARDIZED,
#                      overweight$SE_STANDARDIZED,
#                      percieved_age$SE_STANDARDIZED,
#                      serum_urate_concentration$SE_STANDARDIZED,
#                      telomere_length$SE_STANDARDIZED,
#                      triglyceride$SE_STANDARDIZED,
#                      two_hour_glucose$SE_STANDARDIZED,
#                      type_2_diabetes$SE_STANDARDIZED,
#                      waist_hip_ratio$SE_STANDARDIZED,
#                      atherosclerosis$SE_STANDARDIZED,
#                      leucine$SE_STANDARDIZED,
#                      isoleucine$SE_STANDARDIZED,
#                      valine$SE_STANDARDIZED,
#                      free_cholesterol_in_very_large_hdl$SE_STANDARDIZED,
#                      free_cholesterol_in_large_hdl$SE_STANDARDIZED,
#                      phospholipid_in_very_large_hdl$SE_STANDARDIZED,
#                      phospholipid_in_large_hdl$SE_STANDARDIZED,
#                      triglyceride_in_very_large_vldl$SE_STANDARDIZED,
#                      triglyceride_in_large_vldl$SE_STANDARDIZED), ncol = 42, byrow = F)
# 
# n_matrix = matrix(c(alzheimers$N,
#                     body_mass_index$N,
#                     chronic_inflammation$N,
#                     class_1_obesity$N,
#                     class_2_obesity$N,
#                     class_3_obesity$N,
#                     coronary_artery_disease$N,
#                     extreme_body_mass_index$N,
#                     extreme_waist_hip_ratio$N,
#                     extrinsic_epigenetic_age_acceleration$N,
#                     fasting_glucose$N,
#                     fasting_insulin$N,
#                     fasting_proinsulin$N,
#                     frailty_index$N,
#                     hba1c$N,
#                     healthspan$N,
#                     high_density_lipoprotien$N,
#                     homab$N,
#                     homair$N,
#                     income$N,
#                     intelligence$N,
#                     intrinsic_epigenetic_age_acceleration$N,
#                     lacunar_stroke$N,
#                     low_density_lipoprotien$N,
#                     overweight$N,
#                     percieved_age$N,
#                     serum_urate_concentration$N,
#                     telomere_length$N,
#                     triglyceride$N,
#                     two_hour_glucose$N,
#                     type_2_diabetes$N,
#                     waist_hip_ratio$N,
#                     atherosclerosis$N,
#                     leucine$N,
#                     isoleucine$N,
#                     valine$N,
#                     free_cholesterol_in_very_large_hdl$N,
#                     free_cholesterol_in_large_hdl$N,
#                     phospholipid_in_very_large_hdl$N,
#                     phospholipid_in_large_hdl$N,
#                     triglyceride_in_very_large_vldl$N,
#                     triglyceride_in_large_vldl$N), ncol = 42, byrow = F)
# 




# beta_matrix = matrix(c(homair$BETA_STANDARDIZED,
#                        homab$BETA_STANDARDIZED,
#                        hba1c$BETA_STANDARDIZED,
#                        fasting_glucose$BETA_STANDARDIZED,
#                        fasting_insulin$BETA_STANDARDIZED,
#                        fasting_proinsulin$BETA_STANDARDIZED,
#                        triglyceride$BETA_STANDARDIZED,
#                        serum_urate_concentration$BETA_STANDARDIZED,
#                        waist_hip_ratio$BETA_STANDARDIZED,
#                        body_mass_index$BETA_STANDARDIZED,
#                        chronic_inflammation$BETA_STANDARDIZED,
#                        valine$BETA_STANDARDIZED,
#                        type_2_diabetes$BETA_STANDARDIZED,
#                        leucine$BETA_STANDARDIZED,
#                        isoleucine$BETA_STANDARDIZED,
#                        healthspan$BETA_STANDARDIZED,
#                        atherosclerosis$BETA_STANDARDIZED,
#                        coronary_artery_disease$BETA_STANDARDIZED,
#                        low_density_lipoprotien$BETA_STANDARDIZED,
#                        high_density_lipoprotien$BETA_STANDARDIZED,
#                        lacunar_stroke$BETA_STANDARDIZED,
#                        two_hour_glucose$BETA_STANDARDIZED
#                        ), ncol = length(traits_to_keep), byrow = F)
# 
# se_matrix = matrix(c(homair$SE_STANDARDIZED,
#                      homab$SE_STANDARDIZED,
#                      hba1c$SE_STANDARDIZED,
#                      fasting_glucose$SE_STANDARDIZED,
#                      fasting_insulin$SE_STANDARDIZED,
#                      fasting_proinsulin$SE_STANDARDIZED,
#                      triglyceride$SE_STANDARDIZED,
#                      serum_urate_concentration$SE_STANDARDIZED,
#                      waist_hip_ratio$SE_STANDARDIZED,
#                      body_mass_index$SE_STANDARDIZED,
#                      chronic_inflammation$SE_STANDARDIZED,
#                      valine$SE_STANDARDIZED,
#                      type_2_diabetes$SE_STANDARDIZED,
#                      leucine$SE_STANDARDIZED,
#                      isoleucine$SE_STANDARDIZED,
#                      healthspan$SE_STANDARDIZED,
#                      atherosclerosis$SE_STANDARDIZED,
#                      coronary_artery_disease$SE_STANDARDIZED,
#                      low_density_lipoprotien$SE_STANDARDIZED,
#                      high_density_lipoprotien$SE_STANDARDIZED,
#                      lacunar_stroke$SE_STANDARDIZED,
#                      two_hour_glucose$SE_STANDARDIZED
#                      ), ncol = length(traits_to_keep), byrow = F)

beta_matrix = matrix(c(homair$BETA,
                       homab$BETA,
                       hba1c$BETA,
                       fasting_glucose$BETA,
                       fasting_insulin$BETA,
                       fasting_proinsulin$BETA,
                       triglyceride$BETA,
                       serum_urate_concentration$BETA,
                       waist_hip_ratio$BETA,
                       body_mass_index$BETA,
                       chronic_inflammation$BETA,
                       valine$BETA,
                       type_2_diabetes$BETA,
                       leucine$BETA,
                       isoleucine$BETA,
                       healthspan$BETA,
                       atherosclerosis$BETA,
                       coronary_artery_disease$BETA,
                       low_density_lipoprotien$BETA,
                       high_density_lipoprotien$BETA,
                       lacunar_stroke$BETA,
                       two_hour_glucose$BETA
), ncol = length(traits_to_keep), byrow = F)

se_matrix = matrix(c(homair$SE,
                     homab$SE,
                     hba1c$SE,
                     fasting_glucose$SE,
                     fasting_insulin$SE,
                     fasting_proinsulin$SE,
                     triglyceride$SE,
                     serum_urate_concentration$SE,
                     waist_hip_ratio$SE,
                     body_mass_index$SE,
                     chronic_inflammation$SE,
                     valine$SE,
                     type_2_diabetes$SE,
                     leucine$SE,
                     isoleucine$SE,
                     healthspan$SE,
                     atherosclerosis$SE,
                     coronary_artery_disease$SE,
                     low_density_lipoprotien$SE,
                     high_density_lipoprotien$SE,
                     lacunar_stroke$SE,
                     two_hour_glucose$SE
), ncol = length(traits_to_keep), byrow = F)

n_matrix = matrix(c(homair$N,
                    homab$N,
                    hba1c$N,
                    fasting_glucose$N,
                    fasting_insulin$N,
                    fasting_proinsulin$N,
                    triglyceride$N,
                    serum_urate_concentration$N,
                    waist_hip_ratio$N,
                    body_mass_index$N,
                    chronic_inflammation$N,
                    valine$N,
                    type_2_diabetes$N,
                    leucine$N,
                    isoleucine$N,
                    healthspan$N,
                    atherosclerosis$N,
                    coronary_artery_disease$N,
                    low_density_lipoprotien$N,
                    high_density_lipoprotien$N,
                    lacunar_stroke$N,
                    two_hour_glucose$N
                    ), ncol = length(traits_to_keep), byrow = F)


print(6)


gip1_effect_sizes = beta_matrix %*% eigenvectors[,1]



print(7)

gip1_effect_size_se = ((se_matrix[,13]^2) + ((beta_matrix[,13]^2)/rowSums(n_matrix)))
gip1_effect_size_se = gip1_effect_size_se * gip1_variance
gip1_effect_size_se = gip1_effect_size_se - ((gip1_effect_sizes^2)/rowSums(n_matrix))
gip1_effect_size_se = sqrt(gip1_effect_size_se)


gip1_effect_sizes = gip1_effect_sizes/sqrt(gip1_variance)
gip1_effect_size_se = gip1_effect_size_se/sqrt(gip1_variance)


print(8)

# 
# 
# list_effect_sizes_gip = list()
# for (i in 1:4){
#   gip1_effect_size_se = sqrt(gip1_variance[,i] * ((se_matrix^2) + ((beta_matrix^2)/n_matrix[,1])) - ((gip1_effect_sizes[,i]^2)/n_matrix[,1]))
# 
#   list_effect_sizes_gip1[[i]] = gip1_effect_size_se/sqrt(gip1_variance)
# }





gip1_dataframe = data.frame(gip1_effect_sizes, gip1_effect_size_se)
colnames(gip1_dataframe) = c("beta", "se")
gip1_dataframe$z = gip1_dataframe$beta/gip1_dataframe$se
gip1_dataframe$p_val = pchisq((gip1_dataframe$z^2),1,low=F)#2*pnorm(-abs(gip1_dataframe$z))
gip1_dataframe$neg_log10_p_val = -1 * log(gip1_dataframe$p_val, base = 10)
gip1_dataframe$rsid = common_snps


#print to see rnage of p values
print(max(gip1_dataframe$neg_log10_p_val, na.rm = TRUE))

print(9)

save(gip1_dataframe, eigenvectors, eigenvalues, temp_correlation, file = '~/scratch/principal_component_GWAS/sladek_gip1/gip1_dataframe.RData')


print("done 10")

