# Case-controls analysis 
## Steps


1. Run quack pipeline. Reads were aligned to hg19 reference genome using BWA and SNPs and indels were called according to GATK best practices recommendation.
2. Run pipeline_rabdomyzer.sh. gVCF file were analyzed and annotated using VEP 
3. Genes in cases and controls were selected using Rscript [select_genes.R](https://github.com/Manuelaio/collapse_analysis/blob/main/select_genes.R) with following filters:
      - FILTER =PASS
      - AC =< 1 
      - CADD > 20
      - gnomAD_genomes_AC and gnomAD_exomes_AC <=1 or not present 
      - Consequence: missense, stop gained, stop lost, splice
4. Genes with strong quality difference in cases and controls cohort was excluded  usig Rscript [coverage_analysis.R](https://github.com/Manuelaio/collapse_analysis/blob/main/coverage_analysis.R)
