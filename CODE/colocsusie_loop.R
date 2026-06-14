setwd("/Users/spencypinto/Desktop/COLOC_SUSIE_184")
library(data.table)
library(coloc)
library(susieR)

rm(list = ls())

############################################################
# SETTINGS
############################################################

my_window <- 250000
p_sig_thresh <- 0.001

N_sle <- 25333
s_sle <- 0.27

N_covid <- 1086211
s_covid <- 0.013

loci_file <- "SuSiE_183_loci_with_signals.txt"

sle_dir   <- "gwas_subset_eur"
covid_dir <- "covid_gwas_subset"
ld_dir    <- "LD_matrix_eur_covid_v_sle"

out_file <- "coloc_susie_SLE_EUR_vs_COVID_EUR_results.tsv"

############################################################
# LOAD LOCI
############################################################

loci <- fread(loci_file)

if ("Locus" %in% names(loci)) {
  loci[, lead_snp := Locus]
} else if ("RSID" %in% names(loci)) {
  loci[, lead_snp := RSID]
} else {
  stop("No Locus or RSID column found")
}

all_results <- list()

############################################################
# MAIN LOOP
############################################################

for (i in seq_len(nrow(loci))) {
  
  lead_snp <- loci$lead_snp[i]
  
  cat("\n====================\n")
  cat("Running:", lead_snp, "\n")
  
  sle_file   <- file.path(sle_dir,   paste0("SLE_meta_eur_", lead_snp, ".txt"))
  covid_file <- file.path(covid_dir, paste0("COVID_meta_eur_", lead_snp, ".txt"))
  
  ld_file  <- file.path(ld_dir, paste0(lead_snp, "_LD.unphased.vcor1"))
  var_file <- paste0(ld_file, ".vars")
  
  if (!file.exists(sle_file)) {
    cat("Missing SLE file\n")
    next
  }
  
  if (!file.exists(covid_file)) {
    cat("Missing COVID file\n")
    next
  }
  
  if (!file.exists(ld_file)) {
    cat("Missing LD file\n")
    next
  }
  
  if (!file.exists(var_file)) {
    cat("Missing LD vars file\n")
    next
  }
  
  tryCatch({
    
    ########################################################
    # READ DATA
    ########################################################
    
    sle <- fread(sle_file)
    covid <- fread(covid_file)
    
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
    
    covid2 <- covid[, .(
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
    
    sle2 <- unique(sle2, by = "snp")
    covid2 <- unique(covid2, by = "snp")
    
    ########################################################
    # MERGE COMMON SNPs
    ########################################################
    
    m <- merge(sle2, covid2, by = "snp", suffixes = c(".sle", ".covid"))
    m <- m[chr.sle == chr.covid]
    
    if (nrow(m) < 20) stop("Too few SNPs after merge")
    
    ########################################################
    # REMOVE PALINDROMIC SNPs
    ########################################################
    
    m <- m[!(
      (A1.sle == "A" & A2.sle == "T") |
        (A1.sle == "T" & A2.sle == "A") |
        (A1.sle == "C" & A2.sle == "G") |
        (A1.sle == "G" & A2.sle == "C")
    )]
    
    ########################################################
    # ALIGN ALLELES
    ########################################################
    
    m <- m[
      (A1.sle == A1.covid & A2.sle == A2.covid) |
        (A1.sle == A2.covid & A2.sle == A1.covid)
    ]
    
    if (nrow(m) < 20) stop("Too few SNPs after allele alignment")
    
    m[
      A1.sle == A2.covid & A2.sle == A1.covid,
      `:=`(
        beta.covid = -beta.covid,
        A1.covid = A1.sle,
        A2.covid = A2.sle
      )
    ]
    
    m[, position := position.sle]
    
    ########################################################
    # WINDOW AROUND SLE LEAD SNP
    ########################################################
    
    lead_idx <- which.min(m$p.sle)
    lead_pos <- m$position[lead_idx]
    
    m <- m[
      position >= lead_pos - my_window &
        position <= lead_pos + my_window
    ]
    
    ########################################################
    # KEEP SNPs WITH P < 0.001 IN BOTH TRAITS
    ########################################################
    
    m <- m[p.sle < p_sig_thresh | p.covid < p_sig_thresh]
    
    if (nrow(m) < 15) stop("Too few SNPs after p-value filtering")
    
    ########################################################
    # CLEAN ALLELES + MAF
    ########################################################
    
    m <- m[
      nchar(A1.sle) == 1 & nchar(A2.sle) == 1 &
        A1.sle %in% c("A","C","G","T") &
        A2.sle %in% c("A","C","G","T")
    ]
    
    m[, MAF.sle := pmin(freq.sle, 1 - freq.sle)]
    m[, MAF.covid := pmin(freq.covid, 1 - freq.covid)]
    
    m <- m[
      !is.na(MAF.sle) &
        !is.na(MAF.covid) &
        MAF.sle > 0 &
        MAF.covid > 0
    ]
    
    if (nrow(m) < 10) stop("Too few SNPs after MAF filtering")
    
    ########################################################
    # READ LD MATRIX + VARS FILE
    ########################################################
    
    ld <- as.matrix(fread(ld_file, header = FALSE))
    
    ld_snps <- fread(var_file, header = FALSE)$V1
    
    if (nrow(ld) != length(ld_snps)) {
      stop("LD matrix and .vars length mismatch")
    }
    
    rownames(ld) <- ld_snps
    colnames(ld) <- ld_snps
    
    common <- intersect(m$snp, ld_snps)
    
    if (length(common) < 10) stop("Too few SNPs common with LD")
    
    m <- m[match(common, snp)]
    ld <- ld[common, common]
    
    if (nrow(m) != nrow(ld)) stop("LD and summary stats mismatch")
    
    diag(ld) <- 1
    ld <- (ld + t(ld)) / 2
    
    ########################################################
    # DATASETS
    ########################################################
    
    D1 <- list(
      beta = m$beta.sle,
      varbeta = m$se.sle^2,
      snp = m$snp,
      position = m$position,
      type = "cc",
      s = s_sle,
      N = N_sle,
      MAF = m$MAF.sle,
      LD = ld
    )
    
    D2 <- list(
      beta = m$beta.covid,
      varbeta = m$se.covid^2,
      snp = m$snp,
      position = m$position,
      type = "cc",
      s = s_covid,
      N = N_covid,
      MAF = m$MAF.covid,
      LD = ld
    )
    
    check_dataset(D1)
    check_dataset(D2)
    
    ########################################################
    # RUN SUSIE
    ########################################################
    
    S1 <- runsusie(D1)
    S2 <- runsusie(
      D2,
      estimate_prior_variance = FALSE,
      scaled_prior_variance = 0.15^2
    )
    
    ########################################################
    # RUN COLOC.SUSIE
    ########################################################
    
    res <- coloc.susie(S1, S2)
    
    res_sum <- as.data.table(res$summary)
    
    res_sum[, `:=`(
      locus = lead_snp,
      n_snps = nrow(m),
      SLE_signals = length(S1$sets$cs),
      COVID_signals = length(S2$sets$cs),
      top_sle_snp = m$snp[which.max(S1$pip)],
      top_sle_pip = max(S1$pip),
      top_covid_snp = m$snp[which.max(S2$pip)],
      top_covid_pip = max(S2$pip)
    )]
    
    if ("hit1" %in% names(res_sum) & "hit2" %in% names(res_sum)) {
      res_sum[, ld_between_hits := mapply(function(a, b) {
        if (a %in% rownames(ld) && b %in% colnames(ld)) {
          ld[a, b]
        } else {
          NA_real_
        }
      }, hit1, hit2)]
      
      res_sum[, r2_between_hits := ld_between_hits^2]
    }
    
    all_results[[lead_snp]] <- res_sum
    
    cat("Finished:", lead_snp, "\n")
    
  }, error = function(e) {
    cat("FAILED:", lead_snp, "\n")
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

############################################################
# SAVE RESULTS
############################################################

final_results <- rbindlist(all_results, fill = TRUE)

fwrite(final_results, out_file, sep = "\t")

cat("\nSaved:", out_file, "\n")
cat("Number of loci with results:", length(all_results), "\n")
cat("Number of rows:", nrow(final_results), "\n")