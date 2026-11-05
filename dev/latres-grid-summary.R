# Generates the term x family table of dev/latres-findings.md from the logs of
# dev/latres-grid.R: Rscript dev/latres-grid-summary.R <outdir>
args <- commandArgs(TRUE)
fs <- list.files(args[[1]], "^grid-.*tsv$", full.names = TRUE)
g <- do.call(rbind, lapply(fs, function(f) {
  x <- readLines(f, warn = FALSE)
  m <- strsplit(x[-1L], "\t", fixed = TRUE)
  m <- lapply(m, function(v) c(v, rep("", 10L - length(v))))
  d <- as.data.frame(do.call(rbind, m), stringsAsFactors = FALSE)
  names(d) <- strsplit(x[1L], "\t", fixed = TRUE)[[1L]]
  d
}))
cat("rows:", nrow(g), "\n\n")
sym <- function(i) {
  b <- g$brms[i]
  if (b == "refuses") return("R")
  c0 <- g$compiled[i]
  ev <- g$evaluates[i]
  if (identical(c0, "TRUE") && identical(ev, "TRUE")) return("B+")
  if (identical(c0, "TRUE")) return("B?")
  if (identical(c0, "FALSE")) return("Bx")
  "B"
}
g$cell <- vapply(seq_len(nrow(g)), sym, "")
g$fr <- ifelse(g$frm == "refuses", "r", "f")
g$both <- paste0(g$cell, "/", g$fr)
terms <- c("ar1t", "ar1f", "ar2t", "ar2f", "ma1t", "ma1f", "ma2t", "ma2f",
           "arma1t", "arma1f", "arma2t", "arma2f", "cosy", "unstr", "fcor")
fam <- unique(g[, c("class", "family")])
tab <- matrix("", nrow(fam), length(terms),
              dimnames = list(fam$family, terms))
for (i in seq_len(nrow(g))) tab[g$family[i], g$term[i]] <- g$both[i]
cat("| family |", paste(terms, collapse = " | "), "|\n")
cat("|---|", paste(rep("---", length(terms)), collapse = "|"), "|\n", sep = "")
for (f in rownames(tab)) {
  cat("|", sub(".*", f, f), "|", paste(tab[f, ], collapse = " | "), "|\n")
}
cat("\nbrms's refusals, in its words, with counts:\n")
r <- g[g$brms == "refuses", ]
print(sort(table(substr(r$brms_error, 1, 90)), decreasing = TRUE))
cat("\nbrms builds and frm() refuses:", sum(g$brms == "builds" & g$frm == "refuses"),
    "\nbrms refuses and frm() fits:", sum(g$brms == "refuses" & g$frm != "refuses"),
    "\nboth refuse:", sum(g$brms == "refuses" & g$frm == "refuses"),
    "\nboth build:", sum(g$brms == "builds" & g$frm != "refuses"), "\n")
cat("\nbrms refuses and frm() fits:\n")
print(g[g$brms == "refuses" & g$frm != "refuses", c("term", "family")],
      row.names = FALSE)
cat("\nfrm()'s refusals where brms refuses too, distinct texts:\n")
print(sort(table(substr(g$frm_error[g$brms == "refuses" & g$frm == "refuses"], 1,
                        100)), decreasing = TRUE))
cat("\nfrm()'s refusals where brms builds, distinct texts:\n")
print(sort(table(substr(g$frm_error[g$brms == "builds" & g$frm == "refuses"], 1,
                        100)), decreasing = TRUE))
