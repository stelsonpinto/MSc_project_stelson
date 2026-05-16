library(data.table)
library(coloc)
library(susieR)
rm(list=ls())

#set window around main signal to run SuSie and COLOC on
my_window <- 500000 ##originally was 250000 
#set p-value threshold for filter SNPs on in eQTL and SLE (will remove SNPs ONLY if larger than this in BOTH data)
p_sig_thresh <- 0.001
#---------------------------
# Supply study-level constants
#---------------------------
# Replace these with your actual values
N_eqtl <- 31684        # example only
N_covid <- 1086211
s_covid  <- 0.013         # example only: proportion of SLE cases = total cases/cases+control

#---------------------------
# Read files
#---------------------------
covid  <- fread("covid_gwas_subset/covid_meta_eur_rs9437.txt")
eqtl <- fread("eqtlgen_B38/eqtl_rs9437_ENSG00000104388_B38.txt")

#---------------------------
# Keep and rename columns
#---------------------------
covid2 <- covid[, .(
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

#Check Duplicate SNPs
covid2  <- unique(covid2,  by = "snp")
eqtl2 <- unique(eqtl2, by = "snp")


#---------------------------
# Merge on shared SNPs
#---------------------------
m <- merge(
  covid2, eqtl2,
  by = "snp",
  suffixes = c(".covid", ".eqtl")
)

# Optional: require same chromosome / position if you want a stricter merge
m <- m[chr.covid == chr.eqtl]


#Check Palindromic SNPs (A/T or C/G)
palindromic <- (m$A1.covid == "A" & m$A2.covid == "T") |
  (m$A1.covid == "T" & m$A2.covid == "A") |
  (m$A1.covid == "C" & m$A2.covid == "G") |
  (m$A1.covid == "G" & m$A2.covid == "C")

m <- m[!palindromic]

#---------------------------
# Align alleles
#---------------------------
# Keep only SNPs where alleles match directly or as reversed coding
direct  <- m$A1.covid == m$A1.eqtl & m$A2.covid == m$A2.eqtl
reverse <- m$A1.covid == m$A2.eqtl & m$A2.covid == m$A1.eqtl

m <- m[direct | reverse]

# recompute reverse AFTER filtering
reverse <- m$A1.covid == m$A2.eqtl & m$A2.covid == m$A1.eqtl

# Flip eQTL beta where alleles are reversed
m[reverse, `:=`(beta.eqtl = -beta.eqtl,
                A1.eqtl = A1.covid,
                A2.eqtl = A2.covid)]

# Use one position column
m[, position := position.covid]

#---------------------------
# Restrict to locus window (IMPORTANT)
#---------------------------

# pick lead SNP from eQTL (or SLE)
lead_idx <- which.min(m$p.covid)
lead_pos <- m$position[lead_idx]

# define window (e.g. ±250kb)
window <- my_window

# subset
m <- m[
  position >= lead_pos - window &
    position <= lead_pos + window
]

#---------------------------
# Keep informative SNPs (ADD HERE)
#---------------------------
m <- m[p.covid < p_sig_thresh | p.eqtl < p_sig_thresh]


cat("Number of SNPs after windowing:", nrow(m), "\n")


#---------------------------
# MAF
#---------------------------
# If freq is effect allele frequency rather than MAF, convert:
m[, MAF.covid  := pmin(freq.covid,  1 - freq.covid)]
m[, MAF.eqtl := pmin(freq.eqtl, 1 - freq.eqtl)]

#protects against meta-analysis files with varying N
# protects against varying N
# keep only SNPs with expected sample size
m <- m[N == N_covid]

# sanity check
stopifnot(nrow(m) > 0)
stopifnot(length(unique(m$N)) == 1)


#keep only plain SNPs.
m <- m[
  nchar(A1.covid) == 1 & nchar(A2.covid) == 1 &
    nchar(A1.eqtl) == 1 & nchar(A2.eqtl) == 1 &
    A1.covid %in% c("A","C","G","T") &
    A2.covid %in% c("A","C","G","T") &
    A1.eqtl %in% c("A","C","G","T") &
    A2.eqtl %in% c("A","C","G","T")
]




#---------------------------
# Load LD matrix
#---------------------------
ld <- as.matrix(read.table("LD_matrix_covid/rs9437_ENSG00000104388_LD.unphased.vcor1", header = FALSE))
ld_snps <- scan("LD_matrix_covid/rs9437_ENSG00000104388_LD.unphased.vcor1.vars", what = "")

# Assign SNP names from PLINK output
rownames(ld) <- ld_snps
colnames(ld) <- ld_snps


#---------------------------
# Align LD and m PROPERLY
#---------------------------

# Keep only SNPs present in LD (preserve LD order)
m <- m[snp %in% ld_snps]

# Reorder m to match LD
m <- m[match(ld_snps, snp)]

# Remove SNPs not in m (if any)
keep <- !is.na(m$snp)
m <- m[keep]
ld <- ld[keep, keep]

# Final sanity check
stopifnot(all(m$snp == rownames(ld)))


# remove SNPs with any NA in LD
#bad <- apply(ld, 1, function(x) any(is.na(x)))
#ld <- ld[!bad, !bad]
#m  <- m[!bad]

#---------------------------
# Build coloc/SuSiE input lists
#---------------------------
D1 <- list(
  beta     = m$beta.covid,
  varbeta  = m$se.covid^2,
  snp      = m$snp,
  position = m$position,
  type     = "cc",
  N        = as.integer(N_covid[1]),
  s        = s_covid,
  MAF      = m$MAF.covid,
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


m2 <- m[!is.na(p.covid) & !is.na(p.eqtl)]
min(c(m2$p.covid, m2$p.eqtl))
m2[
  which.min(pmin(p.covid, p.eqtl))
]
#---------------------------
# Check input structure
#---------------------------
check_dataset(D1)
check_dataset(D2)

# This is specifically useful before SuSiE if LD may be misaligned
check_alignment(D1)
check_alignment(D2)


#---------------------------
# Run SuSiE
#---------------------------
S1 <- runsusie(D1)
#S2 <- runsusie(D2)
S2 <- runsusie(D2,estimate_prior_variance = FALSE,scaled_prior_variance = 0.15^2)

#check results for D1
#number of signals
length(S1$sets$cs) #sets returned (wit internal filter)
susie_get_cs(S1) #all credible sets
#get SNP names in sets
cs <- susie_get_cs(S1)
snps <- m$snp[cs$cs$L1]
snps
#get top SNP
idx <- cs$cs$L1
pip <- S1$pip[idx]
best <- idx[which.max(pip)]
cat("Top SNP:", m$snp[best], "\n")
cat("PIP:", max(pip), "\n")
#get all SNPS in credible set
cs <- susie_get_cs(S1)
idx <- cs$cs$L1
cs_snps <- data.table(
  snp = m$snp[idx],
  position = m$position[idx],
  pip = S1$pip[idx]
)
cs_snps[order(-pip)]


#check results for D2
#number of signals
length(S2$sets$cs) #sets returned (wit internal filter)
susie_get_cs(S2) #all credible sets
#get SNP names in sets
cs <- susie_get_cs(S2)
snps <- m$snp[cs$cs$L1]
snps
#get top SNP
idx <- cs$cs$L1
pip <- S2$pip[idx]
best <- idx[which.max(pip)]
cat("Top SNP:", m$snp[best], "\n")
cat("PIP:", max(pip), "\n")
#get all SNPS in crredible set
cs <- susie_get_cs(S2)
idx <- cs$cs$L1
cs_snps <- data.table(
  snp = m$snp[idx],
  position = m$position[idx],
  pip = S2$pip[idx]
)
cs_snps[order(-pip)]

#---------------------------
# Run coloc.susie
#---------------------------
res <- coloc.susie(S1, S2)

# Inspect signal-pair summary
res$summary


print(res$summary)
print(res$summary[order(-res$summary$PP.H4.abf), ])

######################################################################
###define required coverage for causal snps credible set##############
######################################################################
my_coverage <- 0.95
#Find rows with high PP.H4
#choose the row and set below to view all probabilities of biing shared causal (SNP.PP.H4)
# Extract SNPs + posterior for signal pair 4
df <- data.frame(
  snp = res$results$snp,
  pp_h4 = res$results$SNP.PP.H4.row1
)
# Remove missing values (important!)
df <- df[!is.na(df$pp_h4), ]
# Order by decreasing posterior
df <- df[order(df$pp_h4, decreasing = TRUE), ]
# View top SNPs
head(df, 10)
#############################################
#######extratc credible set##################
#############################################
# Compute cumulative sum
df$cum_pp <- cumsum(df$pp_h4)
# Keep SNPs until cumulative posterior reaches 0.95
cs95 <- df[df$cum_pp <= my_coverage, ]
# If you want to include the SNP that crosses 0.95:
cs95 <- df[seq_len(which(df$cum_pp >= my_coverage)[1]), ]
# View credible set
cs95


#########
##plots##
#########
plot(
  m$position,
  -log10(m$p.sle),
  pch = 19,
  xlab = "Position",
  ylab = "-log10(p)",
  main = "SLE locus"
)

plot(
  m$position,
  -log10(m$p.eqtl),
  pch = 19,
  xlab = "Position",
  ylab = "-log10(p)",
  main = "eQTL locus"
)



