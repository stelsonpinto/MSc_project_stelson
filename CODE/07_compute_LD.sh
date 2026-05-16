#!/bin/bash
#SBATCH --job-name=liftover_plink
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/liftover_%j.out


mkdir -p LD_matrix_prep
mkdir -p LD_matrix
Module load plink2


awk 'NR>1 {print $1, $2, $7}' ST3_with_gene_names.tsv | while read rsid chr gene; do

  prefix=${rsid}_${gene}

  if [[ -f LD_ref_data/eqtl_${prefix}_B38_snps.txt && -f LD_ref_data/eqtl_${prefix}_B38_alleles.txt ]]; then

    ./plink2 \
      --pfile chr_data/all_hg38_chr${chr}_EUR_ref_nodup vzs \
      --extract LD_ref_data/eqtl_${prefix}_B38_snps.txt \
      --ref-allele force LD_ref_data/eqtl_${prefix}_B38_alleles.txt 2 1 \
      --make-pgen vzs \
      --out LD_matrix_prep/${prefix}_aligned

    ./plink2 \
      --pfile LD_matrix_prep/${prefix}_aligned vzs \
      --maf 0.02 \
      --r-unphased square ref-based \
      --out LD_matrix/${prefix}_LD

  else
    echo "Skipping ${prefix} (missing files)"
  fi

done
