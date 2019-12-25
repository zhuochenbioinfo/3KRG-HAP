# 3KRG-HAP

Haplotype marker based on 3,000 Rice Genomes for subpopulation inference

***NOTICE: This is a highly specialized pipeline and may be awkward for other uses. But some scripts in this pipeline may be usable for other situations.***


# Pipeline


**1. SNP pruning**

`perl SNP_pruning.r2.pl --in 3K-RG.vcf --out 3K-RG.geno`

This step yielded the 3K-SNP dataset. Named as 3K-SNP.geno in the following steps.


**2. haplotype construction for 3K-RG**

`perl geno_to_binhap.pl --in 3K-RG.geno --out 3K-HAP.haplotype`

This step yielded the 3K-RG haplotype file for each window. Named as 3K-HAP.haplotype


**3. NAF-score calculation in 3K-RG**

Need a tab-delimited list of 3K-RG sample names and subpopulation assignment. Named as 3K-RG.sample_list

`perl haplotype_to_subtype_standard.pl 3K-HAP.haplotype 3K-RG.sample_list 3K-HAP.haplotype.NAF_score`

This step yielded a NAF-scores for each subpopulation on each haplotype. Named as 3K-RG.haplotype.NAF_score


**4. genotyping in custom population**

Make a bed or interval file of SNPs in 3K-SNP dataset as required in GATK UnifiedGenotyper

Perform genotyping using GATK UnifiedGenotyper with these parameters: --L 3K-SNP.bed or --L 3K-SNP.intervals and --output_mode EMIT_ALL_SITES.

This step yielded VCF format genotype file of 3K-SNPs in custom population. Named as custom.vcf.


**5. haplotype construction in custom population**

Make a varlist for haplotype construction with this command:

`cat 3K-SNP.geno | grep -v "^#"| awk '{bin=int(($2-1)/10000);name=sprintf("%05d",bin);print $1"\t"$2"\t"$3"\t"$4"\t"$1"_"name}' > 3K-SNP.varlist`

`perl gatk_vcf_to_haplotype_with_varlist.pl --vcf custom.vcf --var 3K-SNP.varlist --out custom.haplotype`


**6. haplotype matching and NAF-score for custom population**

`perl classify_sample_haplotype_score.pl custom.haplotype 3K-RG.haplotype.NAF_score outpath/`

This step output NAF-score for each 10 kb window of each sample. Named as sample.NAF.

Then merge NAF-score for each 100 kb window.

`perl scan_haplotype_stdratio.pl sample.NAF 10kb_window.bed sample.bin_NAF`


**7. window subpopulation assignment**

`perl dissect_rice_bin.pl sample.bin_NAF > sample.bin_NAF`


**8. plotting sample binmap**
`Rscript draw_bin.rice.R chr.len sample.bin_NAF sample.bin_NAF.pdf`

