#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
Casi_tg=args[1]
Controlli_tg = args[2]
outpath= args[3]

Casi.Topp= read.delim(Casi_tg, header = T, sep="\t",fill = TRUE ,stringsAsFactors = F)
Casi.syn= read.delim("/work/emanuela.iovino/collapse_analysis/new_analysis/analysis_02_04_21/output_toppgene/CPO_SYN_FILTER.txt")
Controlli.Topp= read.delim(Controlli_tg, header = T, sep="\t",fill = TRUE ,stringsAsFactors = F)

MF1=Casi.Topp[Casi.Topp$Category=="GO: Molecular Function",]
MF2=Casi.syn[Casi.syn$Category=="GO: Molecular Function",]
MF3=Controlli.Topp[Controlli.Topp$Category=="GO: Molecular Function",]

MF.casi=subset(MF1, select=c(Name,q.value.FDR.B.H))
MF.casi$logp = -log10(MF.casi[,2])
MF.casi$group='clp.case'

MF.controlli= subset(MF3,  select=c(Name,q.value.FDR.B.H))
MF.controlli$logp = -log10(MF.controlli[,2])
MF.controlli$group='controls'


MF.syn= subset(MF2,  select=c(Name,q.value.FDR.B.H))
MF.syn$logp = -log10(MF.syn[,2])
MF.syn$group='controls.syn'

library(dplyr)
library(ggplot2)


MF.casi.f= MF.casi %>% slice(1:20)
mf.name=as.list(MF.casi.f$Name)
MF.fil.syn=MF.syn %>% dplyr::filter(Name %in% mf.name)
MF.fil.controls=MF.controlli %>% dplyr::filter(Name %in% mf.name)
tot.MF= rbind(MF.casi.f, MF.fil.syn, MF.fil.controls)

tot.MF %>%
  ggplot(aes(x=forcats::fct_reorder(Name,logp),y=logp, fill=group)) +
  geom_col() + 
  coord_flip() 


png(file=paste0(outpath, "/Molecular_Function.png"))
mf= tot.MF %>%
  ggplot(aes(x=forcats::fct_reorder(Name,logp), y=logp, fill=group)) +
  geom_col(position="dodge") + 
  coord_flip() + geom_hline(yintercept = 1, size = 0.5)+
  labs(x="GO: Molecular Function, first 20 GO", y="-log(pvalue)")

mf + scale_fill_manual(values=c("#F8766D", "#00BA38"))
dev.off()
################################################################################


MF1=Casi.Topp[Casi.Topp$Category=="GO: Biological Process",]
MF2=Casi.syn[Casi.syn$Category=="GO: Biological Process",]
MF3=Controlli.Topp[Controlli.Topp$Category=="GO: Biological Process",]

MF.casi=subset(MF1, select=c(Name,q.value.FDR.B.H))
MF.casi$logp = -log10(MF.casi[,2])
MF.casi$group='clp.case'

MF.controlli= subset(MF3,  select=c(Name,q.value.FDR.B.H))
MF.controlli$logp = -log10(MF.controlli[,2])
MF.controlli$group='controls'


MF.syn= subset(MF2,  select=c(Name,q.value.FDR.B.H))
MF.syn$logp = -log10(MF.syn[,2])
MF.syn$group='controls.syn'

MF.casi.f= MF.casi %>% slice(1:25)


mf.name=as.list(MF.casi.f$Name)
MF.fil.syn=MF.syn %>% dplyr::filter(Name %in% mf.name)
MF.fil.controls=MF.controlli %>% dplyr::filter(Name %in% mf.name)
tot.MF= rbind(MF.casi.f, MF.fil.syn, MF.fil.controls)


png(file=paste0(outpath,"/Biological_Process.png"))
tot.MF %>%
  ggplot(aes(x=forcats::fct_reorder(Name,logp), y=logp, fill=group)) +
  geom_col(position="dodge") + 
  coord_flip() + geom_hline(yintercept = 1, size = 0.5)+
  labs(x="GO: Biological Process", y="-log(pvalue)") + ggtitle("Biological Process") + scale_fill_manual(values=c("#F8766D", "#00BA38"))
dev.off()


###############################################################################
#MF1=Casi.Topp[Casi.Topp$Category=="GO: Cellular Component",]
MF1=Casi.Topp[Casi.Topp$Category=="GO: Cellular Component",]
MF2=Casi.syn[Casi.syn$Category=="GO: Cellular Component",]
MF3=Controlli.Topp[Controlli.Topp$Category=="GO: Cellular Component",]

MF.casi=subset(MF1, select=c(Name,q.value.FDR.B.H))
MF.casi$logp = -log10(MF.casi[,2])
MF.casi$group='clp.case'

MF.controlli= subset(MF3,  select=c(Name,q.value.FDR.B.H))
MF.controlli$logp = -log10(MF.controlli[,2])
MF.controlli$group='controls'


MF.syn= subset(MF2,  select=c(Name,q.value.FDR.B.H))
MF.syn$logp = -log10(MF.syn[,2])
MF.syn$group='controls.syn'

MF.casi.f= MF.casi %>% slice(1:25)


mf.name=as.list(MF.casi.f$Name)
MF.fil.syn=MF.syn %>% dplyr::filter(Name %in% mf.name)
MF.fil.controls=MF.controlli %>% dplyr::filter(Name %in% mf.name)
tot.MF= rbind(MF.casi.f, MF.fil.syn, MF.fil.controls)


png(file=paste0(outpath,"/Cellular_component.png"))
tot.MF %>%
  ggplot(aes(x=forcats::fct_reorder(Name,logp), y=logp, fill=group)) +
  geom_col(position="dodge") + 
  coord_flip() + geom_hline(yintercept = 1, size = 0.5)+
  labs(x="GO: Cellular Component", y="-log(pvalue)") + ggtitle("GO: Cellular Component") + scale_fill_manual(values=c("#F8766D", "#00BA38", "#C77CFF"))
dev.off()

################################################################################
MF1=Casi.Topp[Casi.Topp$Category=="Pathway",]
MF2=Casi.syn[Casi.syn$Category=="Pathway",]
MF3=Controlli.Topp[Controlli.Topp$Category=="Pathway",]

MF.casi=subset(MF1, select=c(Name,q.value.FDR.B.H))
MF.casi$logp = -log10(MF.casi[,2])
MF.casi$group='clp.case'

MF.controlli= subset(MF3,  select=c(Name,q.value.FDR.B.H))
MF.controlli$logp = -log10(MF.controlli[,2])
MF.controlli$group='controls'


MF.syn= subset(MF2,  select=c(Name,q.value.FDR.B.H))
MF.syn$logp = -log10(MF.syn[,2])
MF.syn$group='controls.syn'

MF.casi.f= MF.casi %>% slice(1:25)


mf.name=as.list(MF.casi.f$Name)
MF.fil.syn=MF.syn %>% dplyr::filter(Name %in% mf.name)
MF.fil.controls=MF.controlli %>% dplyr::filter(Name %in% mf.name)
tot.MF= rbind(MF.casi.f, MF.fil.syn, MF.fil.controls)

png(file=paste0(outpath,"/patway.png"))

tot.MF %>%
  ggplot(aes(x=forcats::fct_reorder(Name,logp), y=logp, fill=group)) +
  geom_col(position="dodge") + 
  coord_flip() + geom_hline(yintercept = 1, size = 0.5)+
  labs(x="Pathway", y="-log(pvalue)") + ggtitle("Pathway") + scale_fill_manual(values=c("#F8766D", "#00BA38", "#C77CFF"))
dev.off()

