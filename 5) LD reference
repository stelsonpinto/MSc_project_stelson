#!/bin/bash
#SBATCH --job-name=liftover_plink
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/liftover_%j.out

module load plink2

plink2 --pfile all_hg38 vzs --max-alleles 2 --remove remove_related_fixed.txt --keep filtered_samples_EUR_fixed.txt --allow-extra-chr --make-pgen vzs --out all_hg38_EUR

plink2 --pfile all_hg38_EUR vzs --ref-from-fa GRCh38_full_analysis_set_plus_decoy_hla.fa.zst --allow-extra-chr --make-pgen vzs --out all_hg38_EUR_ref

plink2 --pfile all_hg38_EUR_ref vzs --rm-dup force-first --make-pgen vzs --allow-extra-chr --out all_hg38_EUR_ref_nodup

plink2 --pfile all_hg38_EUR_ref_nodup vzs --allow-extra-chr --make-just-pvar --out all_hg38_EUR_ref_nodup_unzipped

