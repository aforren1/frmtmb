# Final check: compat-proposals.tsv against the findings' proposals block,
# row for row, and the nl sweep's proposals against the TSV.
p <- read.delim("dev/equiv-data/compat-proposals.tsv", stringsAsFactors = FALSE)
s <- read.delim("dev/equiv-data/nl-sweep.tsv", stringsAsFactors = FALSE)
fx <- readLines("dev/equiv-findings.md")
b <- grep("^<!-- equiv-proposals begin -->$", fx)
e <- grep("^<!-- equiv-proposals end -->$", fx)
rows <- grep("^\\| ", fx[(b + 1):(e - 1)], value = TRUE)[-1]
cells <- lapply(strsplit(sub("^\\| ", "", sub(" \\|$", "", rows)),
                " | ", fixed = TRUE), trimws)
tab <- do.call(rbind, lapply(cells, function(x) x[1:5]))
cat("TSV rows:", nrow(p), " findings rows:", nrow(tab), "\n")
pair <- gsub("\\|", "|", tab[, 2], fixed = TRUE)
k_md <- paste(tab[, 1], pair, tab[, 3], tab[, 4], tab[, 5])
k_tsv <- paste(ifelse(is.na(p$id), "", p$id),
               paste(p$feature_a, "x", p$feature_b), p$status_now,
               p$proposal, p$verdict)
cat("rows in findings not in TSV:", sum(!k_md %in% k_tsv),
    "; rows in TSV not in findings:", sum(!k_tsv %in% k_md), "\n")
if (any(!k_md %in% k_tsv)) print(k_md[!k_md %in% k_tsv])
if (any(!k_tsv %in% k_md)) print(k_tsv[!k_tsv %in% k_md])
cat("unique pairs:", length(unique(paste(pmin(p$feature_a, p$feature_b),
                                         pmax(p$feature_a, p$feature_b)))),
    "\n")
m <- match(paste("nl x", s$family), paste(p$feature_a, "x", p$feature_b))
d <- !is.na(m) & s$proposal != p$proposal[m]
cat("sweep vs TSV proposal mismatches:", sum(d), "\n")
print(table(p$proposal))
