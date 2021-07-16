#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dir.controls<- args[1]
dir.casi<- args[2]
lov_cvg<- args[3]

#________________load_controls_____________#

ff <- list.files(dir.controls, pattern = ".*coverage.sample_gene_summary", recursive = TRUE, full.names = TRUE)

out_r= do.call(rbind, lapply(ff , function(x)  data.frame(basename(dirname(x)), read.table(x,dec=".", stringsAsFactors = FALSE))))
colnames(out_r)=c("sample", "Genes","total_coverage", "average_coverage", "total_cvg", "mean_cvg","granular_Q1","granular_median", "granular_Q3","cv_above_20")

#list of high-covered genes

out= out_r[out_r$granular_Q3!=">500",]
library(dplyr)
str(out)
out$average_coverage=as.numeric(out$average_coverage)
out$total_coverage=as.numeric(out$total_coverage)
out$total_cvg=as.numeric(out$total_cvg)
out$mean_cvg=as.numeric(out$mean_cvg)
out$granular_Q1=as.numeric(out$granular_Q1)
out$granular_median=as.numeric(out$granular_median)
out$granular_Q3=as.numeric(out$granular_Q3)
out$cv_above_20 =as.numeric(out$cv_above_20)

gene_controlli= out %>%
  dplyr::group_by(Genes) %>%
  dplyr::summarize(mean = mean(average_coverage), md= median(average_coverage),target= mean(cv_above_20))

high_qul= out_r[out_r$granular_Q3==">500",]

#_________________load cases______________________#

ff.casi <- list.files(dir.casi, pattern = ".*coverage.sample_gene_summary", recursive = TRUE, full.names = TRUE)


out_casi= do.call(rbind, lapply(ff.casi , function(x)  data.frame(basename(dirname(x)), read.table(x,dec=".", stringsAsFactors = FALSE))))
colnames(out_casi)=c("sample", "Genes","total_coverage", "average_coverage", "total_cvg", "mean_cvg","granular_Q1","granular_median", "granular_Q3","cv_above_20")
#list of  high-covered genes

outC= out_casi[out_casi$granular_Q3!=">500",]

high_qul_casi= out_casi[out_casi$granular_Q3==">500",]
hig_tot= rbind(high_qul_casi,high_qul)
high= hig_tot[c(2)]

library(dplyr)
str(outC)
outC$average_coverage=as.numeric(outC$average_coverage)
outC$total_coverage=as.numeric(outC$total_coverage)
outC$total_cvg=as.numeric(outC$total_cvg)
outC$mean_cvg=as.numeric(outC$mean_cvg)
outC$granular_Q1=as.numeric(outC$granular_Q1)
outC$granular_median=as.numeric(outC$granular_median)
outC$granular_Q3=as.numeric(outC$granular_Q3)
outC$cv_above_20 =as.numeric(outC$cv_above_20)

gene_casi= outC %>%
  dplyr::group_by(Genes) %>%
  dplyr::summarize(mean_casi = mean(average_coverage), md_casi= median(average_coverage), target_casi= mean(cv_above_20))

val_covergae= dplyr::inner_join(gene_casi, gene_controlli)

val_covergae$ev_mean_target= with(val_covergae, ifelse(val_covergae$target > val_covergae$target_casi, "major", "minor"))


val_covergae$diff <- with(val_covergae, ifelse(ev_mean_target %in% "major", val_covergae$target - val_covergae$target_casi, val_covergae$target_casi - val_covergae$target))


#       Genes to be removed:
#       1 Genes with a mean of cov < 20x in 80% of target 


low_cov_cc=val_covergae[val_covergae$target_casi< 80 |val_covergae$target<80,]
vec_g= low_cov_cc$Genes
geneNAME_lc_= low_cov_cc[c(1)]


#     2 Gens with a target difference > 10 % between cases and controls 

hig_cv= val_covergae %>%
  dplyr::filter(target_casi > 80 | target > 80 )

diff= hig_cv %>% dplyr::filter(!Genes %in% vec_g)
geni_diff_target_cc= diff[diff$diff>15,]
geni_diff_target_cc_= geni_diff_target_cc[c(1)]

#   make a list with low-reliable genes 

list2= geni_diff_target_cc[geni_diff_target_cc$Genes,]
lista_tot= rbind(geneNAME_lc_,geni_diff_target_cc_, high)

write.table(lista_tot, file= lov_cvg, quote=F, col.names = F, row.names = F)





