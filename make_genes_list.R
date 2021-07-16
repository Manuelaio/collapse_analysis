#Info:this script make input file for toppgene and return some summary infomration about sample. The input file of singleton_genes.R will be the output of make_df.R 
#input controls= genes het and hom merged with cat command made by make_df.R


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
controlliSani=args[1]
casi_path = args[2]
low_cv= args[3]
output_controls=args[4]
output_casi= args[5]

library(dplyr)
library(ggplot2)


controlliSani=read.table(controlliSani, stringsAsFactors = F)
colnames(controlliSani)= c("sample_controlli", "CHROM", "POS", "REF",
                           "ALT","Consequence","impact","AC","gnomAD_genomes_AC","gnomAD_exomes_AC", "ExAC_pLI", "SYMBOL","VARIANT_CLASS")

samples=as.data.frame(table(controlliSani$sample))
boxplot(samples$Freq)

singleton= controlliSani %>%
  dplyr::filter(is.na(gnomAD_genomes_AC)) %>%
  dplyr::filter(is.na(gnomAD_exomes_AC)) %>%
  dplyr::filter(VARIANT_CLASS=="SNV") %>%
  dplyr::filter(impact !="LOW")


tab1=as.data.frame(table(singleton$Consequence))
c=prop.table(tab1$Freq)*100
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
c.d<-paste0(specify_decimal(c, 2),"%")
lab<- rep("controls",nrow(tab1))
tab.controlli= cbind(lab,tab1, c.d)
colnames(tab.controlli)= c("group","consequence","Freq","percentage")
library(ggplot2)
# Barplot
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

controlli=ggplot(data = tab.controlli, aes(x = "", y = Freq, fill = consequence)) + 
  geom_bar(stat = "identity") +
  # geom_text(#aes(label = percentage), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") + blank_theme + ggtitle ("singleton variants in controls")
controlli + scale_fill_manual(values=c("#7CAE00","#00BFC4", "#F8766D", "#F8766D"))

gene.controls=data.frame()
for(m in singleton$sample){
  gene.controls=data.frame(unique(singleton$SYMBOL))}



################################
#load file 



casi= read.table(casi_path, stringsAsFactors = F)


colnames(casi)= c("sample", "CHROM", "POS", "REF",
                  "ALT","Consequence","impact","AC","gnomAD_genomes_AC","gnomAD_exomes_AC", "ExAC_pLI", "SYMBOL", "VARIANT_CLASS")

#check how many variants are in samples 

samples=as.data.frame(table(casi$sample))
boxplot(samples$Freq) 


singleton.casi= casi_no_outlier %>%
  dplyr::filter(is.na(gnomAD_genomes_AC)) %>%
  dplyr::filter(is.na(gnomAD_exomes_AC)) %>%
  dplyr::filter(VARIANT_CLASS=="SNV") %>%
  dplyr::filter(impact !="LOW")


tab1=as.data.frame(table(singleton.casi$Consequence))
c=prop.table(tab1$Freq)*100
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
c.d<-paste0(specify_decimal(c, 2),"%")
lab.c<- rep("cases",nrow(tab1))
tab.Case= cbind(lab.c,tab1, c.d)
colnames(tab.Case)= c("group","consequence","Freq","percentage")
library(ggplot2)
# Barplot
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
#COLOR= "#F8766D","#7CAE00","#00BFC4", "#C77CFF"
casi= ggplot(data = tab.Case, aes(x = "", y = Freq, fill = consequence)) + 
  geom_bar(stat = "identity") +
  # geom_text(aes(label = percentage), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") + blank_theme + ggtitle ("singleton variants in casi")
casi + scale_fill_manual(values=c( "#F8766D","#7CAE00","#00BFC4"))

#identificare varianti comuni a casi e controlli
gene.casi=data.frame()
for(m in singleton.casi$sample){
  gene.casi=data.frame(unique(singleton.casi$SYMBOL))}

df= inner_join(singleton, singleton.casi)
geni_da_rimuovere=df$SYMBOL 

#rimovere le varianti comuni a casi e controlli nella gene list dei controlli

gene.controls.filter= gene.controls %>% 
  dplyr::filter(!unique.singleton.SYMBOL. %in% geni_da_rimuovere)



low_cv=read.table(low_cv, stringsAsFactors = F)
vec= low= low_cv$V1

# rimuovere i geni della lista low-coverage dalla lista dei cotnrolli

gene.controls.filter2= gene.controls.filter %>% 
  dplyr::filter(!unique.singleton.SYMBOL. %in% vec)



write.table(gene.controls.filter2, file=output_controls, quote = F,
            row.names = F, col.names = F)


#_casi
#rimovere le varianti comuni a casi e controlli nella gene list dei controlli

gene.casi.filter= gene.casi %>% 
  dplyr::filter(!unique.singleton.casi.SYMBOL. %in% geni_da_rimuovere)

gene.casi.filter2= gene.casi.filter %>% 
  dplyr::filter(!unique.singleton.casi.SYMBOL. %in% vec)



write.table(gene.casi.filter2, file=output_casi, quote = F,
            row.names = F, col.names = F)

