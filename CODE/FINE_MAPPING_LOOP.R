library(data.table)
library(susieR)

rm(list = ls())

#---------------------------
# Parameters
#---------------------------
N_gwas <- 25333
L_max <- 10

#---------------------------
# Output directory
#---------------------------
dir.create("susie_results_eur_sle", showWarnings = FALSE)

#---------------------------
# GWAS files
#---------------------------
gwas_files <- list.files(
  "gwas_subset_eur",
  pattern = "^SLE_meta_eur_.*\\.txt$",
  full.names = TRUE
)

summary_results <- list()

#---------------------------
# Loop over loci
#---------------------------
for(i in seq_along(gwas_files)) {
  
  gwas_file <- gwas_files[i]
  
  locus <- sub(
    "^SLE_meta_eur_(.*)\\.txt$",
    "\\1",
    basename(gwas_file)
  )
  
  cat("\n====================================\n")
  cat("Processing:", locus, "\n")
  
  ld_file <- file.path(
    "LD_matrix_eur_sle",
    paste0(locus, "_LD.unphased.vcor1")
  )
  
  vars_file <- file.path(
    "LD_matrix_eur_sle",
    paste0(locus, "_LD.unphased.vcor1.vars")
  )
  
  if (!file.exists(ld_file) || !file.exists(vars_file)) {
    cat("Missing LD files for", locus, "\n")
    next
  }
  
  tryCatch({
    
    #---------------------------
    # Read data
    #---------------------------
    sumstats <- fread(gwas_file)
    
    R <- as.matrix(read.table(ld_file))
    
    ld_snps <- fread(vars_file, header = FALSE)
    colnames(ld_snps) <- "rsid"
    
    #---------------------------
    # Match SNP order
    #---------------------------
    idx <- match(ld_snps$rsid, sumstats$RSID)
    
    if (sum(is.na(idx)) > 0) {
      cat(
        "Skipping", locus,
        "-", sum(is.na(idx)),
        "LD SNPs absent from GWAS\n"
      )
      next
    }
    
    sumstats <- sumstats[idx, ]
    
    stopifnot(
      all(sumstats$RSID == ld_snps$rsid)
    )
    
    #---------------------------
    # Calculate Z scores
    #---------------------------
    z <- sumstats$b_eur / sumstats$se_eur
    
    if (length(z) != nrow(R)) {
      cat(
        "Skipping", locus,
        "- dimension mismatch\n"
      )
      next
    }
    
    #---------------------------
    # Ensure LD matrix symmetric
    #---------------------------
    R <- (R + t(R)) / 2
    
    #---------------------------
    # Run SuSiE
    #---------------------------
    fit <- susie_rss(
      z = z,
      R = R,
      n = N_gwas,
      L = L_max
    )
    
    #---------------------------
    # SNP-level results
    #---------------------------
    pip <- data.frame(
      RSID = sumstats$RSID,
      Chr  = sumstats$Chr,
      Pos  = sumstats$Pos,
      PIP  = fit$pip,
      CS   = NA_character_
    )
    
    #credible set
    if (!is.null(fit$sets$cs)) {
      
      for (j in seq_along(fit$sets$cs)) {
        
        pip$CS[fit$sets$cs[[j]]] <- paste0("CS", j)
        
      }
      
    }
    
    # Sort by PIP
    pip <- pip[order(-pip$PIP), ]
    
    #---------------------------
    # Save per-locus results
    #---------------------------
    fwrite(
      pip,
      file.path(
        "susie_results_eur_sle",
        paste0(locus, "_susie_pip.txt")
      ),
      sep = "\t"
    )
    
    saveRDS(
      fit,
      file.path(
        "susie_results_eur_sle",
        paste0(locus, "_susie_fit.rds")
      )
    )
    
    #---------------------------
    # Summary statistics
    #---------------------------
    top_snp <- pip$RSID[1]
    top_pip <- pip$PIP[1]
    top_chr <- pip$Chr[1]
    top_pos <- pip$Pos[1]
    
    n_signals <- ifelse(
      is.null(fit$sets$cs),
      0,
      length(fit$sets$cs)
    )
    
    cs95_size <- ifelse(
      is.null(fit$sets$cs),
      0,
      length(unique(unlist(fit$sets$cs)))
    )
    
    summary_results[[length(summary_results) + 1]] <- data.frame(
      Index = i,
      Locus = locus,
      Chr = top_chr,
      Pos = top_pos,
      Top_SNP = top_snp,
      Top_PIP = top_pip,
      Signals = n_signals,
      CS95_Size = cs95_size,
      stringsAsFactors = FALSE
    )
    
    cat(
      "Finished:",
      locus,
      "| Top SNP:", top_snp,
      "| PIP:", round(top_pip, 4),
      "| Signals:", n_signals,
      "| CS95 size:", cs95_size,
      "\n"
    )
    
  }, error = function(e) {
    
    cat(
      "ERROR:",
      locus,
      "-",
      e$message,
      "\n"
    )
    
  })
  
}

#---------------------------
# Save summary file
#---------------------------
summary_df <- rbindlist(summary_results)

fwrite(
  summary_df,
  "SuSiE_183_loci_summary.txt",
  sep = "\t"
)

cat(
  "\nCompleted",
  nrow(summary_df),
  "loci\n"
)

print(head(summary_df))
