library(stats)
omim = read.table("list.omim.txt")
hpo = read.table("list.human_pheno_ontology.txt")
genotology = read.table("list.genontology.txt")
genecard = read.table("list.geneCard.txt")
malacard = read.table("list.malacard.txt")

#convert in list 
omim_list = as.character(unlist(omim))
hpo_list =as.character(unlist(hpo))
genotology_list = as.character(unlist(genotology))
genecard_list = as.character(unlist(genecard))

malacard_list = as.character(unlist(malacard))


# intersected genes 

intesection_genes = Reduce(intersect,list(hpo_list,omim_list, malacard_list, genotology_list,genecard_list)) #25 genes


###select enriched genes in case 

genes_fisher= read.table("../new_analysis/cpo_analysis/all_sample/collapse_all/collapse_fisher.txt", header=T)
enrch = data.frame(genes_fisher[genes_fisher$p_val< 0.05 & genes_fisher$n_varianti_controlli<=10,]$GENES)
enrch_list = as.character(unlist(enrch))

genes = list(intesection_genes,enrch_list)

#number of genes enriched in list of cpo genes 

observed <-length(Reduce(intersect,list(intesection_genes,enrch_list)))
observed

#Hypergeometric test 

hyper_pval <- phyper(
  q=observed-1, m=length(genes[[1]]),
  n=nrow(genes_fisher)- length(genes[[1]]),
  k=length(genes[[2]]), lower.tail=FALSE
)

hyper_pval # 7.81e-03



