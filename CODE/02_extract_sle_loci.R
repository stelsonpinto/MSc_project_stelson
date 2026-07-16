###create slurm script

#!/bin/bash
#SBATCH --job-name=locus_extract
#SBATCH --output=gwas_subset/locus_extract_%j.out
#SBATCH --error=gwas_subset/locus_extract_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

# Initialize conda
# Activate environment
conda activate r_env
# Go to working directory
# Create output directory
mkdir -p gwas_subset

# Run R script
Rscript extracting_gwas.R


###R script

rm(list=ls())

library(data.table)
library(ggplot2)

############################
# DEFINE VARIABLES
############################
my_study <- "eur"   # "all", "eur", "eas"

if(my_study == "eur"){my_sample_size <- 25333}
if(my_study == "eas"){my_sample_size <- 267034}
if(my_study == "all"){my_sample_size <- 296197}

in_file <- paste("metal_",my_study,"_finalized_rsid.assoc",sep="")

############################
# LOAD GWAS SUMMARY STATS
############################
GWAS_summary_stats <- read.delim(in_file, sep = "\t", stringsAsFactors = FALSE)

GWAS_summary_stats$Chr <- toupper(as.character(GWAS_summary_stats$Chr))

# Convert chromosome X labels to 23
GWAS_summary_stats$Chr[
  GWAS_summary_stats$Chr %in% c("X","CHRX","CHR23")
] <- "23"

if(!"Pos" %in% colnames(GWAS_summary_stats)){
  if("BP" %in% colnames(GWAS_summary_stats)){
    names(GWAS_summary_stats)[names(GWAS_summary_stats)=="BP"] <- "Pos"
  } else {
    stop("No position column found (Pos or BP)")
  }
}

GWAS_summary_stats <- GWAS_summary_stats[GWAS_summary_stats$MarkerName!="",]
GWAS_summary_stats <- GWAS_summary_stats[GWAS_summary_stats$RSID!="",]

names(GWAS_summary_stats)[names(GWAS_summary_stats)=="MarkerName"] <- "SNP"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="Allele1"] <- "A1"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="Allele2"] <- "A2"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="Freq1"] <- "freq"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="Effect"] <- "b"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="StdErr"] <- "se"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="Pval"] <- "p"

GWAS_summary_stats$N <- ifelse(
  is.na(GWAS_summary_stats$TotalSampleSize),
  my_sample_size,
  GWAS_summary_stats$TotalSampleSize
)

############################
# LOAD ST3 REFERENCE FILE 
############################
st3 <- read.delim("ST3_EUR_summary_stats.tsv", stringsAsFactors = FALSE)

st3 <- st3[, c("RSID","chr","start","end")]
st3$chr <- as.character(st3$chr)

st3$chr <- toupper(st3$chr)

# Convert chromosome X labels to 23
st3$chr[
  st3$chr %in% c("X","CHRX","CHR23")
] <- "23"

st3 <- unique(st3)


############################
# LOOP OVER LOCI
############################
for(i in seq_len(nrow(st3))){

  rsid <- st3$RSID[i]
  mychr <- st3$chr[i]
  my_lower <- st3$start[i]
  my_upper <- st3$end[i]

  cat("Processing:", rsid, "\n")

  locus_index <- GWAS_summary_stats$Chr == mychr &
                 GWAS_summary_stats$Pos >= my_lower &
                 GWAS_summary_stats$Pos <= my_upper 

  if(sum(locus_index) == 0){
    message("No variants in subset: ", rsid)
    next
  }

  SLE_meta_for_database <- data.frame(
    SNP = GWAS_summary_stats$SNP[locus_index],
    A1 = GWAS_summary_stats$A1[locus_index],
    A2 = GWAS_summary_stats$A2[locus_index],
    freq = GWAS_summary_stats$freq[locus_index],
    b = GWAS_summary_stats$b[locus_index],
    se = GWAS_summary_stats$se[locus_index],
    p = GWAS_summary_stats$p[locus_index],
    N = GWAS_summary_stats$N[locus_index],
    Chr = GWAS_summary_stats$Chr[locus_index],
    Pos = GWAS_summary_stats$Pos[locus_index],
    RSID = GWAS_summary_stats$RSID[locus_index]
  )

  colnames(SLE_meta_for_database)[4] <- paste("freq_", my_study, sep="")
  colnames(SLE_meta_for_database)[5] <- paste("b_", my_study, sep="")
  colnames(SLE_meta_for_database)[6] <- paste("se_", my_study, sep="")
  colnames(SLE_meta_for_database)[7] <- paste("p_", my_study, sep="")
  colnames(SLE_meta_for_database)[8] <- paste("N_", my_study, sep="")

  ############################
  # SAVE FILE (RSID-BASED)
  ############################
  out_file <- paste0(
    "gwas_subset/subset_data/SLE_meta_",
    my_study, "_", rsid, ".txt"
  )

  write.table(SLE_meta_for_database, file=out_file,
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


cat("All loci processed.\n")
