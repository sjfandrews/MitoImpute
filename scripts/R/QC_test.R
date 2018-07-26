require(rmarkdown)

MtPlatforms = "BDCHP-1X10-HUMANHAP240S_11216501_A-b37"
s = "~/GitCode/MitoImpute/scripts/R/MT_imputation_QC.Rmd"
wgs_map = "~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map"
wgs_ped = "~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped"
wgs_vcf = "~/GitCode/MitoImpute/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz"
typ_map = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, ".map")
typ_ped = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, ".ped")
typ_vcf = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, ".vcf.gz")
imp_map = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, "_imputed.map")
imp_ped = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, "_imputed.ped")
imp_vcf = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, "_imputed.vcf")
imp_info = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, "_imputed_info")

output = paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms, "/chrMT_1kg_", MtPlatforms, "_mtImputed_QC.html")

wd = "~/GitCode/MitoImpute/"

output_dir= paste0("~/GitCode/MitoImpute/DerivedData/ThousandGenomes/", MtPlatforms)
info_cut = "0"

render(s,
       output_file = output,
       output_dir = output_dir,
       params = list(rwd = wd,
                     info_cut = info_cut,
                     wgs_map = wgs_map,
                     wgs_ped = wgs_ped,
                     wgs_vcf = wgs_vcf,
                     typ_map = typ_map,
                     typ_ped = typ_ped,
                     typ_vcf = typ_vcf,
                     imp_map = imp_map,
                     imp_ped = imp_ped,
                     imp_vcf = imp_vcf,
                     imp_info = imp_info))

#rmarkdown::render(rmarkdown::render(${s}, output_file = ${output}, output_dir = ${output_dir}, params = list(rwd = ${rwd}, info.cut = ${info_cut}, wgs.map = ${wgs_map}, wgs.ped = ${wgs_ped}, wgs.vcf = ${wgs_vcf}, typ.map = {input.typ_map}, typ.ped = ${typ_ped}, typ.vcf = ${typ_vcf}, imp.map = ${imp_map}, imp.ped = ${imp_ped}, imp.vcf= ${imp_vcf}, imp.info = ${input.imp_info}))

