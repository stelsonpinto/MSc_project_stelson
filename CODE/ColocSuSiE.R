#R script to run coloc.susie on GWAS and eqtl data
#Run though to lines 168 -> "################create LD ref data" 
#Then create LD ref data using plink and return to line 189 "# Load LD matrix"
getwd()
setwd("/Users/spencypinto/Desktop/COLOC_SUSIE_184")
library(data.table)
library(coloc)
library(susieR)
rm(list=ls())

#set window around main signal to run SuSie and COLOC on
my_window <- 250000 ##originally was 250000 
#set p-value threshold for filter SNPs on in eQTL and SLE (will remove SNPs ONLY if larger than this in BOTH data)
p_sig_thresh <- 0.001
#---------------------------
# Supply study-level constants
#---------------------------
# Replace these with your actual values
N_sle <- 25333
s_sle <- 0.27

# COVID
N_covid <- 1086211
s_covid <- 0.013
# example only: proportion of SLE cases

#---------------------------
# Read files
#---------------------------
sle  <- fread("gwas_subset_eur/SLE_meta_eur_rs2298428.txt")
covid <- fread("covid_gwas_subset/COVID_meta_eur_rs2298428.txt")

#---------------------------
# Keep and rename columns
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
  p        = p_eur,
  N        = N_sle
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
  p        = p_eur,
  N        = N_covid
)]

#Check Duplicate SNPs
sle2  <- unique(sle2,  by = "snp")
covid2 <- unique(covid2, by = "snp")


#---------------------------
# Merge on shared SNPs
#---------------------------
m <- merge(
  sle2, covid2,
  by = "snp",
  suffixes = c(".sle", ".covid")
)

# Optional: require same chromosome / position if you want a stricter merge
m <- m[chr.sle == chr.covid]


#Check Palindromic SNPs (A/T or C/G)
palindromic <- (m$A1.sle == "A" & m$A2.sle == "T") |
  (m$A1.sle == "T" & m$A2.sle == "A") |
  (m$A1.sle == "C" & m$A2.sle == "G") |
  (m$A1.sle == "G" & m$A2.sle == "C")

m <- m[!palindromic]

#---------------------------
# Align alleles
#---------------------------
# Keep only SNPs where alleles match directly or as reversed coding
direct  <- m$A1.sle == m$A1.covid & m$A2.sle == m$A2.covid
reverse <- m$A1.sle == m$A2.covid & m$A2.sle == m$A1.covid

m <- m[direct | reverse]

# recompute reverse AFTER filtering
reverse <- m$A1.sle == m$A2.covid & m$A2.sle == m$A1.covid

# Flip eQTL beta where alleles are reversed
m[reverse, `:=`(beta.covid = -beta.covid,
                A1.covid = A1.sle,
                A2.covid = A2.sle)]

# Use one position column
m[, position := position.sle]

#---------------------------
# Restrict to locus window (IMPORTANT)
#---------------------------

# pick lead SNP from eQTL (or SLE)
lead_idx <- which.min(m$p.sle)
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
m <- m[p.sle < p_sig_thresh | p.covid < p_sig_thresh]


cat("Number of SNPs after windowing:", nrow(m), "\n")


#---------------------------
# MAF
#---------------------------
# If freq is effect allele frequency rather than MAF, convert:
m[, MAF.sle  := pmin(freq.sle,  1 - freq.sle)]
m[, MAF.covid := pmin(freq.covid, 1 - freq.covid)]

#protects against meta-analysis files with varying N
# protects against varying N
# keep only SNPs with expected sample size
m <- m[
  N.covid == N_covid &
    N.sle == N_sle
]

cat("After N filter:", nrow(m), "\n")

# sanity check
stopifnot(nrow(m) > 0)



#keep only plain SNPs.
m <- m[
  nchar(A1.sle) == 1 & nchar(A2.sle) == 1 &
    nchar(A1.covid) == 1 & nchar(A2.covid) == 1 &
    A1.sle %in% c("A","C","G","T") &
    A2.sle %in% c("A","C","G","T") &
    A1.covid %in% c("A","C","G","T") &
    A2.covid %in% c("A","C","G","T")
]


#---------------------------
# Load LD matrix
#---------------------------
ld <- as.matrix(read.table("LD_matrix_eur_covid_v_sle/rs2298428_LD.unphased.vcor1", header = FALSE))
ld_snps <- scan("LD_matrix_eur_covid_v_sle/rs2298428_LD.unphased.vcor1.vars", what = "")

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
  beta     = m$beta.sle,
  varbeta  = m$se.sle^2,
  snp      = m$snp,
  position = m$position,
  type     = "cc",
  N        = as.integer(N_sle[1]),
  s        = s_sle,
  MAF      = m$MAF.sle,
  LD       = ld
)

D2 <- list(
  beta     = m$beta.covid,
  varbeta  = m$se.covid^2,
  snp      = m$snp,
  position = m$position,
  type     = "quant",
  N        = N_covid,
  MAF      = m$MAF.covid,
  LD       = ld
)


m2 <- m[!is.na(p.sle) & !is.na(p.covid)]
min(c(m2$p.sle, m2$p.covid))
m2[
  which.min(pmin(p.sle, p.covid))
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


