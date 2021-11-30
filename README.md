# Case-controls analysis 
## Steps


1. Run quack pipeline. Reads were aligned to hg19 reference genome using BWA and SNPs and indels were called according to GATK best practices recommendation.
2. Run pipeline_rabdomyzer.sh. gVCF file were analyzed and annotated using VEP 
3. PCA analysis performed with [pca_plot.R](https://github.com/Manuelaio/collapse_analysis/blob/main/pca_plot.R) and remove outside samples
4. Genes in cases and controls were selected using Rscript [select_genes.R](https://github.com/Manuelaio/collapse_analysis/blob/main/select_genes.R) with following filters:
      - FILTER =PASS
      - AC =< 1 
      - CADD > 20
      - gnomAD_genomes_AC and gnomAD_exomes_AC <=1 or not present 
      - Consequence: missense, stop gained, stop lost, splice
      
4. Genes with strong quality difference in cases and controls cohort was excluded  usig Rscript [coverage_analysis.R](https://github.com/Manuelaio/collapse_analysis/blob/main/coverage_analysis.R)
5. Two suitables genes lists were prepared with [make_genes_list.R](https://github.com/Manuelaio/collapse_analysis/blob/main/make_genes_list.R) and uploaded on [Toppgene](https://toppgene.cchmc.org) for enrichment analysis 
6. Analysis of cases and controls enrichment analysis was made with [ToppFun.R](https://github.com/Manuelaio/collapse_analysis/blob/main/ToppFun.R)

![alt text](https://github.com/Manuelaio/collapse_analysis/blob/main/GEA.png)

