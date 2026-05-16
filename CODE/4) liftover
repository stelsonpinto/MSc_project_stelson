library(GenomicRanges)
library(rtracklayer)

chain <- import.chain("hg19ToHg38.over.chain")

files <- list.files(
  "eqtlgen",
  pattern = "^eqtl_.*\\.txt$",
  full.names = TRUE
)

for (infile in files) {

  outfile <- sub("\\.txt$", "_B38.txt", infile)

  cat("Processing:", infile, "\n")

  dat <- read.table(infile, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

  dat$Chr <- as.character(dat$Chr)
  dat$Chr[dat$Chr == "23"] <- "X"
  dat$Chr[dat$Chr == "24"] <- "Y"
  dat$Chr[dat$Chr == "25"] <- "XY"
  dat$Chr[dat$Chr == "26"] <- "M"

  gr0 <- GRanges(
    seqnames = paste0("chr", dat$Chr),
    ranges   = IRanges(start = dat$BP, width = 1),
    strand   = "*"
  )
  names(gr0) <- seq_len(length(gr0))

  lo_list <- liftOver(gr0, chain)

  keep     <- which(lengths(lo_list) == 1)
  gr1      <- unlist(lo_list[keep])
  orig_idx <- as.integer(names(gr1))

  chr_out <- as.character(seqnames(gr1))
  chr_out <- sub("^chr", "", chr_out)

  out <- data.frame(
    SNP         = dat$SNP[orig_idx],
    Chr         = chr_out,
    BP          = start(gr1),
    A1          = dat$A1[orig_idx],
    A2          = dat$A2[orig_idx],
    Freq        = dat$Freq[orig_idx],
    Probe       = dat$Probe[orig_idx],
    Probe_Chr   = dat$Probe_Chr[orig_idx],
    Probe_bp    = dat$Probe_bp[orig_idx],
    Gene        = dat$Gene[orig_idx],
    Orientation = dat$Orientation[orig_idx],
    b           = dat$b[orig_idx],
    SE          = dat$SE[orig_idx],
    p           = dat$p[orig_idx],
    stringsAsFactors = FALSE
  )

  write.table(out, file = outfile, sep = "\t",
              quote = FALSE, row.names = FALSE)

  cat("Done:", infile,
      "| kept:", nrow(out),
      "| dropped:", nrow(dat) - nrow(out), "\n\n")
}
