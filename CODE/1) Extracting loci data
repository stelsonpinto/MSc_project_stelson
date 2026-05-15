library(data.table)
library(readxl)
install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

#---------------------------
# 1. Load ST3
#---------------------------
st3 <- as.data.table(
  read_excel("LOCI_LIST.xlsx", sheet = "ST3", skip = 1)
)

#---------------------------
# 2. Select relevant columns
#---------------------------
loci <- st3[, .(
  RSID = RSID,
  chr  = CHR,
  pos  = POS38,
  start = START,
  end   = STOP,
  nearby_gene = `NEARBY GENE(S)`
)]

#---------------------------
# 3. Split multiple genes into rows
#---------------------------
loci <- loci[, .(
  gene_symbol = unlist(strsplit(nearby_gene, ","))
), by = .(RSID, chr, pos, start, end)]

# trim whitespace
loci[, gene_symbol := trimws(gene_symbol)]

# remove empty entries
loci <- loci[gene_symbol != "" & !is.na(gene_symbol)]

#---------------------------
# 4. Map gene symbols → Ensembl IDs
#---------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = unique(loci$gene_symbol),
  mart = mart
)

#---------------------------
# 5. Merge mapping
#---------------------------
loci <- merge(
  loci,
  mapping,
  by.x = "gene_symbol",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

#---------------------------
# 6. Clean final table
#---------------------------
loci_final <- loci[, .(
  RSID,
  chr,
  pos,
  start,
  end,
  gene_symbol,
  gene_id = ensembl_gene_id
)]

# remove rows without mapping
loci_final <- loci_final[!is.na(gene_id)]

#---------------------------
# 7. Save
#---------------------------
fwrite(loci_final, "ST3_with_genes.tsv", sep = "\t")

cat("Done. Rows:", nrow(loci_final), "\n")
