############################
# DEFINE VARIABLES
############################
my_study <- "eas"

# Total GWAS sample size
# Cases = 13,769 (eur) 794 (eas)
# Controls = 1,072,442 (eur) (4862 eas)
# Total N = 1086211 (eur) 5656 (eas)
my_sample_size <- 5656

in_file <- "COVID19_HGI_A2_ALL_eas_leave23andme_20220403.tsv"

############################
# LOAD GWAS SUMMARY STATS
############################
GWAS_summary_stats <- fread(in_file, data.table = FALSE)

############################
# STANDARDIZE COLUMN NAMES
############################

# Chromosome + position
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="#CHR"] <- "Chr"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="POS"]   <- "Pos"

# SNP identifiers
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="rsid"] <- "RSID"

# Alleles
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="ALT"] <- "A1"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="REF"] <- "A2"

# Frequency
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="all_meta_AF"] <- "freq"

# Effect sizes
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="all_inv_var_meta_beta"]   <- "b"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="all_inv_var_meta_sebeta"] <- "se"
names(GWAS_summary_stats)[names(GWAS_summary_stats)=="all_inv_var_meta_p"]      <- "p"

############################
# ADD TOTAL SAMPLE SIZE
############################
GWAS_summary_stats$N <- my_sample_size

############################
# CLEAN / FORMAT DATA
############################

GWAS_summary_stats$Chr <- gsub(
  "^CHR",
  "",
  toupper(as.character(GWAS_summary_stats$Chr))
)

GWAS_summary_stats$Pos <- as.numeric(GWAS_summary_stats$Pos)

# Remove missing SNPs
GWAS_summary_stats <- GWAS_summary_stats[
  GWAS_summary_stats$SNP != "" &
    GWAS_summary_stats$RSID != "",
]

# Remove duplicated variants
GWAS_summary_stats <- GWAS_summary_stats[
  !duplicated(
    paste(
      GWAS_summary_stats$Chr,
      GWAS_summary_stats$Pos,
      GWAS_summary_stats$A1,
      GWAS_summary_stats$A2
    )
  ),
]

############################
# LOAD ST3 REFERENCE FILE
############################
st3 <- fread(
  "COLOC_SUSIE_184/ST3_EUR_summary_stats.tsv",
  data.table = FALSE
)

st3 <- st3[, c("RSID","chr","start","end")]

st3$chr <- gsub(
  "^CHR",
  "",
  toupper(as.character(st3$chr))
)

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
  
  locus_index <- !is.na(GWAS_summary_stats$Chr) &
    !is.na(GWAS_summary_stats$Pos) &
    GWAS_summary_stats$Chr == mychr &
    GWAS_summary_stats$Pos >= my_lower &
    GWAS_summary_stats$Pos <= my_upper
  
  if(sum(locus_index, na.rm = TRUE) == 0){
  
    message("No variants in subset: ", rsid)
    next
  }
  
  COVID_meta_for_database <- data.frame(
    SNP  = GWAS_summary_stats$SNP[locus_index],
    A1   = GWAS_summary_stats$A1[locus_index],
    A2   = GWAS_summary_stats$A2[locus_index],
    freq = GWAS_summary_stats$freq[locus_index],
    b    = GWAS_summary_stats$b[locus_index],
    se   = GWAS_summary_stats$se[locus_index],
    p    = GWAS_summary_stats$p[locus_index],
    N    = GWAS_summary_stats$N[locus_index],
    Chr  = GWAS_summary_stats$Chr[locus_index],
    Pos  = GWAS_summary_stats$Pos[locus_index],
    RSID = GWAS_summary_stats$RSID[locus_index]
  )
  
  ############################
  # RENAME STUDY-SPECIFIC COLUMNS
  ############################
  colnames(COVID_meta_for_database)[4] <- paste0("freq_", my_study)
  colnames(COVID_meta_for_database)[5] <- paste0("b_", my_study)
  colnames(COVID_meta_for_database)[6] <- paste0("se_", my_study)
  colnames(COVID_meta_for_database)[7] <- paste0("p_", my_study)
  colnames(COVID_meta_for_database)[8] <- paste0("N_", my_study)
  
  ############################
  # SAVE FILE
  ############################
  out_file <- paste0(
    "covid_gwas_subset_eas/COVID_meta_",
    my_study,
    "_",
    rsid,
    ".txt"
  )
  
  write.table(
    COVID_meta_for_database,
    file = out_file,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  ############################
  # PLOT
  ############################
  p_col <- paste0("p_", my_study)
  
  df <- COVID_meta_for_database
  
  df[[p_col]] <- as.numeric(df[[p_col]])
  
  df[[p_col]][
    is.na(df[[p_col]]) | df[[p_col]] <= 0
  ] <- .Machine$double.xmin
  
  df$logp <- -log10(df[[p_col]])
  df$Pos  <- as.numeric(df$Pos)
  
  top_i <- which.max(df$logp)
  top_df <- df[top_i, , drop = FALSE]
  
  mean_N <- round(
    mean(df[[paste0("N_", my_study)]], na.rm = TRUE),
    0
  )
  
  my_title <- paste0(
    "Chr ",
    mychr,
    "  ",
    rsid,
    "  ",
    my_study,
    " (N≈",
    format(mean_N, big.mark=","),
    ")"
  )
  
  gg <- ggplot(df, aes(x = Pos, y = logp)) +
    geom_point(size = 0.9, alpha = 0.7) +
    geom_hline(
      yintercept = -log10(5e-8),
      linetype = "dashed"
    ) +
    labs(
      title = my_title,
      x = "Position (bp)",
      y = expression(-log[10](p))
    ) +
    theme_bw()
  
  gg <- gg +
    geom_point(
      data = top_df,
      size = 1.4
    ) +
    geom_text(
      data = top_df,
      aes(label = RSID),
      vjust = -0.7,
      size = 3
    )
  
  print(gg)
}

cat("All loci processed.\n")
