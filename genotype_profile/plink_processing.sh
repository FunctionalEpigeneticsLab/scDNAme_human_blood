#!/usr/bin/bash

# prepare reference
vcf=1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.vcf.gz

plink --vcf $vcf --make-bed --out ref

cut -f2 ref.bim | sort > ref.snps

# check missingness
plink --bfile ref --missing --out ref

# prepare data for the samples
vcf=sciMET_BC_all_10_donors_merged.unfiltered.vcf.gz

plink --vcf $vcf --double-id --snps-only just-acgt --make-bed --out target

plink --bfile target --set-missing-var-ids 'chr@_#_$1_$2' --make-bed --out target.fixed

plink --bfile target.fixed --snps-only --autosome \
      --geno 0 \
      --make-bed \
      --out target.filtered

cut -f2 target.filtered.bim | sort > target.snps

# keep only intersecting SNPs
comm -12 ref.snps target.snps > common.snps

# extract common SNPs
plink --bfile target.fixed --extract common.snps --make-bed --out target.common
plink --bfile ref --extract common.snps --make-bed --out ref.common

# merge reference and target
echo target.common > merge_list.txt
plink --bfile ref.common --merge-list merge_list.txt --make-bed --out merged

# basic SNP filtering (for PCA only)
plink --bfile merged --snps-only --autosome --maf 0.05 --geno 0.01 --hwe 1e-6 --make-bed --out merged.qc

# LD pruning
plink \
  --bfile merged.qc \
  --indep-pairwise 200 50 0.1 \
  --out merged.prune

# extract pruned SNPs:
plink \
  --bfile merged.qc \
  --extract merged.prune.prune.in \
  --make-bed \
  --out merged.pruned
