#!/usr/bin/env Rscript

#Info: this script make as output a dataframe with filtered variants het and hom given a path of rabdomyzer. The filters are: AC<=1, CADD >20, FILTER=PASS, gnomad AC genome and exome >=1. 
# For other filtration use singleton_genes.R

#Run: Rscript select_genes.R /path/rabdomazer /path_output/genes.het.txt   /path_output/genes.hom.txt



args = commandArgs(trailingOnly=TRUE)
dir_rab<- args[1]
het_output= args[2]
hom_output= args[3]

ff <- list.files(dir_rab, pattern = ".*.xlsx", recursive = TRUE, full.names = TRUE)
out=lapply(ff, function(x)  data.frame(basename(dirname(x)), readxl::read_xlsx(x, sheet = "het_calls")))
library(dplyr)

#filter_function 

apply.filters= function(list_df){
  list_df%>%
    dplyr::filter(Consequence== "missense_variant" | Consequence=="stop_gained" | Consequence=="stop_lost" |  Consequence=="splice_region_variant") %>%
    dplyr::filter(AC <=1) %>%
    dplyr::filter(CADD_PHRED >=20 | CADD_PHRED=="." ) %>%
    dplyr::filter(FILTER == "PASS") %>%
    #dplyr::filter(ExAC_pLI >=0.9) %>%
    dplyr::filter(gnomAD_genomes_AC <=1 | is.na(gnomAD_genomes_AC)) %>%
    dplyr::filter(gnomAD_exomes_AC <=1 | is.na(gnomAD_exomes_AC))
}

df.filter= lapply(out, apply.filters)  


for(i in 1:length(df.filter)){
  assign(paste("df", i, sep = ""), df.filter[[i]])
}

all.genes= list()

for(i in 1:length(df.filter)){
  df= df.filter[[i]] %>%
    dplyr::select(basename.dirname.x.., X.CHROM, POS, REF, ALT,Consequence,IMPACT,AC,gnomAD_genomes_AC,gnomAD_exomes_AC, ExAC_pLI, SYMBOL, VARIANT_CLASS)
  all.genes= rlist::list.append(all.genes,df)
}

GENES<- plyr::ldply(all.genes, data.frame)
genes_for_ont.het= unique(GENES)
write.table(genes_for_ont.het, file=het_output, quote= F, row.names = F, col.names = F)

####

out.hom=lapply(ff , function(x)  data.frame(basename(dirname(x)), readxl::read_xlsx(x, sheet = "hom_calls")))
#library(dplyr)

#Applicazione dei filtri 

apply.filters.hom= function(list_df){
  list_df%>%
    dplyr::filter(Consequence== "missense_variant" | Consequence=="stop_gained" | Consequence=="stop_lost" | Consequence=="splice_region_variant") %>%
    dplyr::filter(AC <=2) %>%
    dplyr::filter(CADD_PHRED >=20) %>%
    dplyr::filter(FILTER == "PASS") %>%
    #dplyr::filter(ExAC_pLI >=0.9) %>%
    dplyr::filter(gnomAD_genomes_AC <=1 | is.na(gnomAD_genomes_AC)) %>%
    dplyr::filter(gnomAD_exomes_AC <=1 | is.na(gnomAD_exomes_AC))
}

df.filter.hom= lapply(out.hom, apply.filters.hom)  


for(i in 1:length(df.filter.hom)){
  assign(paste("df", i, sep = ""), df.filter.hom[[i]])
}

homo.genes= list()

for(i in 1:length(df.filter.hom)){
  df= df.filter.hom[[i]] %>%
    dplyr::select(basename.dirname.x.., X.CHROM, POS, REF, ALT,Consequence,IMPACT,AC,gnomAD_genomes_AC,gnomAD_exomes_AC, ExAC_pLI, SYMBOL, VARIANT_CLASS)
  homo.genes= rlist::list.append(homo.genes,df)
}

GENES.hom<- plyr::ldply(homo.genes, data.frame)
genes_for_ont.hom= unique(GENES.hom)
write.table(genes_for_ont.hom, file=hom_output, quote= F, row.names = F, col.names = F)
