#!/usr/bin/env Rscript
# gridss_to_typed.R
#
# Converts a GRIDSS breakend (BND) VCF into simple typed SVs (DEL/DUP/INV/INS),
# following the standard GRIDSS "simple-event-annotation" recipe using the
# StructuralVariantAnnotation Bioconductor package.
#
# Same-chromosome breakpoint pairs are classified by orientation + size:
#   - DEL : deletion-like (++ orientation, ref span removed)
#   - DUP : tandem-duplication-like
#   - INV : inversion-like
#   - INS : insertion (large inserted sequence at the breakend)
# Inter-chromosomal breakpoints (translocations) are dropped — the Q100 DEL/INS
# truth does not contain them.
#
# Usage:
#   Rscript gridss_to_typed.R <input_gridss.vcf.gz> <output_typed.vcf>

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(StructuralVariantAnnotation)
})

args <- commandArgs(trailingOnly = TRUE)
in_vcf  <- args[1]
out_vcf <- args[2]

vcf <- readVcf(in_vcf, "hg38")

# Keep only PASS breakends (GRIDSS quality filter).
vcf <- vcf[rowRanges(vcf)$FILTER == "PASS"]

# Extract breakpoint pairs (gr has one row per breakend; partners linked).
gr <- breakpointRanges(vcf)

# Classify simple event type from breakend orientation + position.
# svtype() from StructuralVariantAnnotation returns DEL/DUP/INV/INS/BND/CTX.
svtype <- simpleEventType(gr)
gr$svtype <- svtype

# Keep only same-chromosome simple types; drop inter-chromosomal (CTX/BND).
keep <- svtype %in% c("DEL", "DUP", "INV", "INS")
gr <- gr[keep]

# Build a minimal typed VCF: one record per breakpoint (use the lower-coord
# partner of each pair to avoid duplicates).
# partner() gives the mate; keep the breakend that starts the event.
is_lower <- start(gr) <= start(partner(gr))
gr1 <- gr[is_lower]

# SVLEN: distance between paired breakends (approx event size).
mate <- partner(gr)[is_lower]
svlen <- abs(start(mate) - start(gr1))
# For DEL svlen is negative by VCF convention.
sign <- ifelse(gr1$svtype == "DEL", -1, 1)

# Write a simple VCF.
con <- file(out_vcf, "w")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##source=GRIDSS-simple-event-annotation",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
  "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
  "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
  "##ALT=<ID=DEL,Description=\"Deletion\">",
  "##ALT=<ID=DUP,Description=\"Duplication\">",
  "##ALT=<ID=INV,Description=\"Inversion\">",
  "##ALT=<ID=INS,Description=\"Insertion\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
), con)

chrom <- as.character(seqnames(gr1))
pos   <- start(gr1)
endp  <- start(mate)
ids   <- names(gr1)
svt   <- gr1$svtype

for (i in seq_along(gr1)) {
  e <- max(pos[i], endp[i])
  line <- sprintf("%s\t%d\t%s\tN\t<%s>\t.\tPASS\tSVTYPE=%s;SVLEN=%d;END=%d",
                  chrom[i], pos[i],
                  ifelse(is.na(ids[i]), paste0("gridss_", i), ids[i]),
                  svt[i], svt[i], as.integer(sign[i] * svlen[i]), e)
  writeLines(line, con)
}
close(con)

cat(sprintf("Wrote %d typed SVs to %s\n", length(gr1), out_vcf))
cat("Type breakdown:\n")
print(table(svt))
