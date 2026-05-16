library(data.table)
rm(list = ls())

#---------------------------
# Parameters
#---------------------------
my_window <- 250000 
p_sig_thresh <- 0.001

N_eqtl <- 31684
N_sle  <- 25333

#---------------------------
# FILES
#---------------------------
eqtl_files <- list.files("eqtlgen_B38", pattern = "eqtl_.*_B38\\.txt$", full.names = TRUE)
gwas_files <- list.files("gwas_subset", 
                         pattern = "SLE_meta_eur_.*\\.txt$", 
                         full.names = TRUE)

st3 <- fread("ST3_EUR_summary_stats.tsv")

dir.create("LD_ref_data", showWarnings = FALSE)

eqtl_files <- eqtl_files[!grepl("eqtl_NA", eqtl_files)]

cat("Number of eQTL files:", length(eqtl_files), "\n")
cat("Number of GWAS files:", length(gwas_files), "\n")

if (length(eqtl_files) == 0) stop("No eQTL files found")
if (length(gwas_files) == 0) stop("No GWAS files found")

#---------------------------
# MAIN LOOP
#---------------------------
for (eqtl_file in eqtl_files) {

  cat("\nProcessing:", eqtl_file, "\n")

  eqtl <- fread(eqtl_file)

  #---------------------------
  # Extract rsID
  #---------------------------
  rsid <- sub("^eqtl_(rs[0-9]+)_.*$", "\\1", basename(eqtl_file))
  cat("Extracted rsid:", rsid, "\n")

  #---------------------------
  # Match GWAS (fixed)
  #---------------------------
  gwas_match <- gwas_files[
    grepl(paste0("_", rsid, "\\.txt$"), gwas_files)
  ]

  if (length(gwas_match) == 0) {
    cat("No GWAS match for", rsid, "\n")
    next
  }

  if (length(gwas_match) > 1) {
    cat("Multiple GWAS matches for", rsid, "\n")
    next
  }

  gwas_file <- gwas_match[1]
  cat("Matched GWAS:", basename(gwas_file), "\n")

  sle <- fread(gwas_file)

  #---------------------------
  # Get locus
  #---------------------------
  locus <- st3[RSID == rsid]

  if (nrow(locus) == 0) {
    cat("No locus found in ST3 for", rsid, "\n")
    next
  }

  #---------------------------
  # FORMAT
  #---------------------------
  sle2 <- sle[, .(
    snp      = RSID,
    chr      = Chr,
    position = Pos,
    A1       = A1,
    A2       = A2,
    freq     = freq_eur,
    beta     = b_eur,
    se       = se_eur,
    p        = p_eur
  )]

  eqtl2 <- eqtl[, .(
    snp      = SNP,
    chr      = Chr,
    position = BP,
    A1       = A1,
    A2       = A2,
    freq     = Freq,
    beta     = b,
    se       = SE,
    p        = p
  )]

  sle2  <- unique(sle2,  by = "snp")
  eqtl2 <- unique(eqtl2, by = "snp")

  #---------------------------
  # REGION FILTER
  #---------------------------
  sle2  <- sle2[position >= locus$start & position <= locus$end]
  eqtl2 <- eqtl2[position >= locus$start & position <= locus$end]

  if (nrow(sle2) < 50) {
    cat("Too few GWAS SNPs in region:", rsid, "\n")
    next
  }

  #---------------------------
  # MERGE
  #---------------------------
  m <- merge(sle2, eqtl2, by = "snp", suffixes = c(".sle", ".eqtl"))
  m <- m[chr.sle == chr.eqtl]

  if (nrow(m) < 20) {
    cat("Too few SNPs after merge:", rsid, "\n")
    next
  }

  #---------------------------
  # REMOVE PALINDROMIC
  #---------------------------
  m <- m[!(
    (A1.sle == "A" & A2.sle == "T") |
    (A1.sle == "T" & A2.sle == "A") |
    (A1.sle == "C" & A2.sle == "G") |
    (A1.sle == "G" & A2.sle == "C")
  )]

  #---------------------------
  # ALIGN (fixed)
  #---------------------------
  m <- m[
    (A1.sle == A1.eqtl & A2.sle == A2.eqtl) |
    (A1.sle == A2.eqtl & A2.sle == A1.eqtl)
  ]

  if (nrow(m) < 20) {
    cat("Too few SNPs after alignment:", rsid, "\n")
    next
  }

  m[
    A1.sle == A2.eqtl & A2.sle == A1.eqtl,
    `:=`(
      beta.eqtl = -beta.eqtl,
      A1.eqtl = A1.sle,
      A2.eqtl = A2.sle
    )
  ]

  m[, position := position.sle]

  #---------------------------
  # WINDOW
  #---------------------------
  lead_idx <- which.min(m$p.sle)
  lead_pos <- m$position[lead_idx]

  m <- m[
    position >= lead_pos - my_window &
    position <= lead_pos + my_window
  ]

  #---------------------------
  # SIGNIFICANT FILTER
  #---------------------------
  m <- m[
    (!is.na(p.sle)  & p.sle  < p_sig_thresh) |
    (!is.na(p.eqtl) & p.eqtl < p_sig_thresh)
  ]

  cat("Number of SNPs after windowing:", nrow(m), "\n")

  if (nrow(m) < 10) {
    cat("Too few SNPs after filtering:", rsid, "\n")
    next
  }

  #---------------------------
  # CLEAN SNPs
  #---------------------------
  m <- m[
    nchar(A1.sle) == 1 & nchar(A2.sle) == 1 &
    nchar(A1.eqtl) == 1 & nchar(A2.eqtl) == 1 &
    A1.sle %in% c("A","C","G","T") &
    A2.sle %in% c("A","C","G","T") &
    A1.eqtl %in% c("A","C","G","T") &
    A2.eqtl %in% c("A","C","G","T")
  ]

  #---------------------------
  # MAF
  #---------------------------
  m[, MAF.sle  := pmin(freq.sle,  1 - freq.sle)]
  m[, MAF.eqtl := pmin(freq.eqtl, 1 - freq.eqtl)]

  m[, N := N_sle]

  #---------------------------
  # LD FILES
  #---------------------------
  region_snps <- unique(m$snp)

  if (length(region_snps) == 0) {
    cat("No SNPs for LD:", rsid, "\n")
    next
  }

  prefix <- paste0(
    "LD_ref_data/",
    sub("\\.txt$", "", basename(eqtl_file))
  )

  write.table(region_snps,
              paste0(prefix, "_snps.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)

  allele_dt <- m[match(region_snps, snp), .(A1.sle, snp)]
  setnames(allele_dt, "A1.sle", "A1")

  write.table(allele_dt,
              paste0(prefix, "_alleles.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)

  cat("LD files written for", basename(prefix), "\n")
}
