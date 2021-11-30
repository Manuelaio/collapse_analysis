
#################################################################################################
#____plot_pca_perferomed_with_EIGENSTRAT#
#################################################################################################


file1= read.table("plink.fam")
file1$new=file1$V1

#sample_list file contain two columns: sample name and groups 

samples= read.table("./samples.list", stringsAsFactors = F, col.names = c("new", "group"))

df= dplyr::inner_join(file1,samples)
fam_plink= df[c(7,7,3,4,5,8)]
write.table(fam_plink, file= "./case_controls.fam", sep= "\t", quote= F, row.names = F, col.names = F)

#file pedsnp
map= read.table("plink.map") 

# zgrep -v "##" recode.cpo.snps-hapmap.vcf.recode.vcf.gz | cut -f 3,4,5 > snp.txt
#snp file contanins snps, REF and ALT

snp= read.table("snp.txt", header= T, stringsAsFactors = F)

snp[snp$ALT=="*",]
colnames(map)[2]= "ID"
#map.ref.alt= merge(map, snp, all = F)
#d=merge(map,unique(snp$ID))

mp= dplyr::left_join(map,unique(snp))

map.ref.alt.uniq=unique(mp)
d <- dplyr::mutate(map.ref.alt.uniq, new = V4/100000000)
map.final= d[,c(1,2,7,4,5,6)]

write.table(map.final, file= "./case_controls.pedsnp", sep= "\t", quote= F, row.names = F, col.names = F)


#__________________________________# out of R 

#smartpca -p pca.par > smartpca.log


#__________________________________#

fn = "./case_controls.pca.evec"
evecDat = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5",
                                     "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))
plot(evecDat$PC1, evecDat$PC2, xlab="PC1", ylab="PC2")
plot(evecDat$PC2, evecDat$PC3, xlab="PC1", ylab="PC2")

d = evecDat[evecDat$Pop=="cpo",]
points(d$PC1, d$PC2, col="red", pch=20)

d = evecDat[evecDat$Pop=="cpo_iran",]
points(d$PC1, d$PC2, col="blue", pch=20)

d = evecDat[evecDat$Pop=="controls_healt",]
points(d$PC1, d$PC2, col="magenta", pch=20)

popgroup <- list(
  CASI= "cpo",
  CASI_IRAN= "cpo_iran",
  CONTROLLI= "controls_healt")

colors <- sapply(levels(evecDat$Pop), function(x) {
  for (i in 1:length(popgroup)) {
    if (x %in% popgroup[[i]])
      return(names(popgroup)[i])
  }
  NA
})
colors <- as.factor(colors)
legend.text <- sapply(levels(colors), function(x) paste(levels(evecDat$Pop)[colors==x], collapse=","))
legend.text

plot(evecDat$PC1, evecDat$PC2, pch=20, cex=0.75, main="Casi-Controlli",
     xlab="eigenvector 1", ylab="eigenvector 2", col=colors[evecDat$Pop])
legend("bottomleft", legend=legend.text, col=1:length(legend.text), pch=19, cex=0.75)
