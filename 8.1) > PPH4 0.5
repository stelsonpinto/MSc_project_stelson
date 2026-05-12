to get everty locus with > PPH4 0.5


library(data.table)
library(coloc)
library(susieR)

rm(list = ls())

######################################################################
######################## SETTINGS #####################################
######################################################################

# Window around lead SNP
my_window <- 250000

# SNP filter threshold
p_sig_thresh <- 0.001

# Study constants
N_eqtl <- 31684
N_sle  <- 25333
s_sle  <- 0.27

######################################################################
######################## LOAD LOCI ####################################
######################################################################

loci <- fread("ST3_with_gene_names.tsv")

# Rename columns
setnames(
  loci,
  old = c("RSID", "gene_id"),
  new = c("lead_snp", "gene")
)

# Remove missing SNPs
loci <- loci[!is.na(lead_snp)]

######################################################################
######################## RESULTS CONTAINER #############################
######################################################################

all_results <- list()

######################################################################
######################## MAIN LOOP ####################################
######################################################################

for(i in 1:nrow(loci)){

  lead_snp <- loci$lead_snp[i]
  gene     <- loci$gene[i]

  cat("\n====================================================\n")
  cat("Running:", lead_snp, gene, "\n")
  cat("====================================================\n")

  out <- tryCatch({

    ############################################################
    # FILE PATHS
    ############################################################

    gwas_file <- paste0(
      "gwas_subset/SLE_meta_eur_",
      lead_snp,
      ".txt"
    )

    eqtl_file <- paste0(
      "eqtlgen_B38/eqtl_",
      lead_snp,
      "_",
      gene,
      "_B38.txt"
    )

    ld_file <- paste0(
      "LD_matrix/",
      lead_snp,
      "_",
      gene,
      "_LD.unphased.vcor1"
    )

    ld_var_file <- paste0(
      "LD_matrix/",
      lead_snp,
      "_",
      gene,
      "_LD.unphased.vcor1.vars"
    )

    ############################################################
    # CHECK FILES EXIST
    ############################################################

    if(!file.exists(gwas_file)){
      stop(paste("Missing GWAS file:", gwas_file))
    }

    if(!file.exists(eqtl_file)){
      stop(paste("Missing eQTL file:", eqtl_file))
    }

    if(!file.exists(ld_file)){
      stop(paste("Missing LD file:", ld_file))
    }

    if(!file.exists(ld_var_file)){
      stop(paste("Missing LD SNP file:", ld_var_file))
    }

    ############################################################
    # READ FILES
    ############################################################

    sle  <- fread(gwas_file)
    eqtl <- fread(eqtl_file)

    ############################################################
    # KEEP / RENAME COLUMNS
    ############################################################

    sle2 <- sle[, .(
      snp      = RSID,
      chr      = Chr,
      position = Pos,
      A1       = A1,
      A2       = A2,
      freq     = freq_eur,
      beta     = b_eur,
      se       = se_eur,
      p        = p_eur,
      N        = N_eur
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

    ############################################################
    # REMOVE DUPLICATES
    ############################################################

    sle2  <- unique(sle2,  by = "snp")
    eqtl2 <- unique(eqtl2, by = "snp")

    ############################################################
    # MERGE
    ############################################################

    m <- merge(
      sle2,
      eqtl2,
      by = "snp",
      suffixes = c(".sle", ".eqtl")
    )

    cat("After merge:", nrow(m), "\n")

    ############################################################
    # SAME CHROMOSOME
    ############################################################

    m <- m[chr.sle == chr.eqtl]

    cat("After chr filter:", nrow(m), "\n")

    ############################################################
    # REMOVE PALINDROMIC SNPs
    ############################################################

    palindromic <- (
      (m$A1.sle == "A" & m$A2.sle == "T") |
      (m$A1.sle == "T" & m$A2.sle == "A") |
      (m$A1.sle == "C" & m$A2.sle == "G") |
      (m$A1.sle == "G" & m$A2.sle == "C")
    )

    m <- m[!palindromic]

    cat("After palindrome removal:", nrow(m), "\n")

    ############################################################
    # ALIGN ALLELES
    ############################################################

    direct <- (
      m$A1.sle == m$A1.eqtl &
      m$A2.sle == m$A2.eqtl
    )

    reverse <- (
      m$A1.sle == m$A2.eqtl &
      m$A2.sle == m$A1.eqtl
    )

    m <- m[direct | reverse]

    reverse <- (
      m$A1.sle == m$A2.eqtl &
      m$A2.sle == m$A1.eqtl
    )

    m[reverse, `:=`(
      beta.eqtl = -beta.eqtl,
      A1.eqtl   = A1.sle,
      A2.eqtl   = A2.sle
    )]

    cat("After allele alignment:", nrow(m), "\n")

    ############################################################
    # POSITION
    ############################################################

    m[, position := position.sle]

    ############################################################
    # WINDOW
    ############################################################

    lead_idx <- which.min(m$p.sle)
    lead_pos <- m$position[lead_idx]

    m <- m[
      position >= lead_pos - my_window &
      position <= lead_pos + my_window
    ]

    cat("After window:", nrow(m), "\n")

    ############################################################
    # P FILTER
    ############################################################

    m <- m[
      p.sle < p_sig_thresh |
      p.eqtl < p_sig_thresh
    ]

    cat("After p filter:", nrow(m), "\n")

    ############################################################
    # MAF
    ############################################################

    m[, MAF.sle  := pmin(freq.sle, 1 - freq.sle)]
    m[, MAF.eqtl := pmin(freq.eqtl, 1 - freq.eqtl)]

    ############################################################
    # SAMPLE SIZE FILTER
    ############################################################

    m <- m[N == N_sle]

    cat("After N filter:", nrow(m), "\n")

    ############################################################
    # KEEP STANDARD SNPs
    ############################################################

    m <- m[
      nchar(A1.sle) == 1 &
      nchar(A2.sle) == 1 &
      nchar(A1.eqtl) == 1 &
      nchar(A2.eqtl) == 1 &
      A1.sle %in% c("A","C","G","T") &
      A2.sle %in% c("A","C","G","T") &
      A1.eqtl %in% c("A","C","G","T") &
      A2.eqtl %in% c("A","C","G","T")
    ]

    cat("After SNP QC:", nrow(m), "\n")

    ############################################################
    # LOAD LD
    ############################################################

    ld <- as.matrix(
      read.table(
        ld_file,
        header = FALSE
      )
    )

    ld_snps <- scan(
      ld_var_file,
      what = ""
    )

    rownames(ld) <- ld_snps
    colnames(ld) <- ld_snps

    ############################################################
    # ALIGN LD
    ############################################################

    m <- m[snp %in% ld_snps]

    cat("After LD overlap:", nrow(m), "\n")

    m <- m[match(ld_snps, snp)]

    keep <- !is.na(m$snp)

    m  <- m[keep]
    ld <- ld[keep, keep]

    cat("After LD alignment:", nrow(m), "\n")

    ############################################################
    # FINAL CHECKS
    ############################################################

    if(nrow(m) < 10){
      stop("Too few SNPs after filtering")
    }

    stopifnot(all(m$snp == rownames(ld)))

    ############################################################
    # BUILD DATASETS
    ############################################################

    D1 <- list(
      beta     = m$beta.sle,
      varbeta  = m$se.sle^2,
      snp      = m$snp,
      position = m$position,
      type     = "cc",
      N        = as.integer(N_sle),
      s        = s_sle,
      MAF      = m$MAF.sle,
      LD       = ld
    )

    D2 <- list(
      beta     = m$beta.eqtl,
      varbeta  = m$se.eqtl^2,
      snp      = m$snp,
      position = m$position,
      type     = "quant",
      N        = N_eqtl,
      MAF      = m$MAF.eqtl,
      LD       = ld
    )

    ############################################################
    # CHECK DATASETS
    ############################################################

    check_dataset(D1)
    check_dataset(D2)

    check_alignment(D1)
    check_alignment(D2)

    ############################################################
    # RUN SUSIE
    ############################################################

    S1 <- runsusie(D1)

    S2 <- runsusie(
      D2,
      estimate_prior_variance = FALSE,
      scaled_prior_variance = 0.15^2
    )

    ############################################################
    # RUN COLOC
    ############################################################

    res <- coloc.susie(S1, S2)

    print(res$summary)

    ############################################################
    # EXTRACT BEST H4
    ############################################################

    if(
      is.null(res$summary) ||
      nrow(res$summary) == 0 ||
      !"PP.H4.abf" %in% colnames(res$summary)
    ){

      best_h4 <- NA

    } else {

      h4_values <- res$summary$PP.H4.abf

      h4_values <- h4_values[
        is.finite(h4_values)
      ]

      if(length(h4_values) == 0){

        best_h4 <- NA

      } else {

        best_h4 <- max(h4_values)

      }

    }

    cat("Best H4 =", best_h4, "\n")

    ############################################################
    # SAVE RESULT
    ############################################################

    result_dt <- data.table(
      locus     = paste0(lead_snp, "_", gene),
      lead_snp  = lead_snp,
      gene      = gene,
      nsnps     = nrow(m),
      PPH4      = best_h4,
      status    = "SUCCESS",
      error     = NA_character_
    )

    result_dt

  }, error = function(e){

    cat("FAILED:", lead_snp, gene, "\n")
    cat("Reason:", e$message, "\n")

    result_dt <- data.table(
      locus     = paste0(lead_snp, "_", gene),
      lead_snp  = lead_snp,
      gene      = gene,
      nsnps     = NA,
      PPH4      = NA,
      status    = "FAILED",
      error     = e$message
    )

    result_dt

  })

  ############################################################
  # STORE RESULT
  ############################################################

  all_results[[length(all_results) + 1]] <- out

}

######################################################################
######################## FINAL OUTPUT #################################
######################################################################

final_results <- rbindlist(
  all_results,
  fill = TRUE
)

############################################################
# SAVE ALL RESULTS
############################################################

fwrite(
  final_results,
  "all_coloc_results.tsv",
  sep = "\t"
)

############################################################
# SAVE HITS
############################################################

hits <- final_results[
  !is.na(PPH4) &
  PPH4 > 0.5
]

fwrite(
  hits,
  "coloc_hits_PPH4_gt_0.5.tsv",
  sep = "\t"
)

######################################################################
######################## SUMMARY ######################################
######################################################################

cat("\nDONE\n")

cat(
  "Total loci tested:",
  nrow(final_results),
  "\n"
)

cat(
  "Successful runs:",
  sum(final_results$status == "SUCCESS"),
  "\n"
)

cat(
  "Failed runs:",
  sum(final_results$status == "FAILED"),
  "\n"
)

cat(
  "Hits with PPH4 > 0.5:",
  nrow(hits),
  "\n"
)
