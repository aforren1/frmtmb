# Lane silent2, item 5207: summary of dev/silent2-ordmeasure.R over the
# controls. Usage: Rscript dev/silent2-ordsum.R <arm> (reads
# dev/silent2-out/ord-<set>-<arm>.tsv for the four sets)
arm <- commandArgs(TRUE)[1]
sets <- c("forty", "sweep", "steep", "wide")
X <- do.call(rbind, lapply(sets, function(s) {
  x <- utils::read.delim(sprintf("dev/silent2-out/ord-%s-%s.tsv", s, arm))
  x$set <- s
  x[, setdiff(names(x), "warn")]
}))
num <- function(v) suppressWarnings(as.numeric(v))
cols <- function(p) X[, grep(paste0("^", p, "[123]$"), names(X)), drop = FALSE]
mx <- function(M) apply(M, 1, function(r) {
  r <- num(r)
  if (all(is.na(r))) NA else max(r, na.rm = TRUE)
})
X$se_t <- mx(cols("se_t"))
X$se_t_nan <- apply(cols("se_t"), 1, function(r) any(is.nan(num(r))))
X$se_b <- num(X$se_b)
cat("fits", nrow(X), "by set:", paste(names(table(X$set)), table(X$set),
                                      collapse = ", "), "\n")
cat("degenerate warning on 0.82.0:", sum(X$degen), "\n")
clean <- !X$degen & X$code == 0
cat("passed (no degenerate warning, code 0):", sum(clean), "\n")
cat("\nlargest threshold SE on passed fits, top 15:\n")
o <- order(-X$se_t[clean])
print(X[clean, ][o[1:15], c("case", "se_t", "se_b", "near1", "near2",
                            "near3", "cover1", "cover2")], row.names = FALSE)
cat("\npassed fits with a NaN threshold SE:\n")
print(X[clean & X$se_t_nan, c("case", "se_t", "se_b", "near1", "near2")],
      row.names = FALSE)
for (cut in c(50, 100, 500, 5e3)) {
  hit <- clean & (X$se_t > cut) %in% TRUE
  cat(sprintf("\nSE cut %g: passed fits flagged %d of %d: %s\n", cut,
              sum(hit), sum(clean), paste(X$case[hit], collapse = " ")))
}
cat("\nlargest coefficient SE on passed fits not flagged at 100, top 8:\n")
k <- clean & !(X$se_t > 100) %in% TRUE
o <- order(-X$se_b[k])
print(X[k, ][o[1:8], c("case", "se_t", "se_b")], row.names = FALSE)
s <- sort(X$se_t[clean & !is.na(X$se_t)])
j <- which.max(diff(log(s)))
cat("\nthe widest gap in the sorted threshold SEs of passed fits:", s[j],
    "->", s[j + 1L], "(ratio", signif(s[j + 1L] / s[j], 3), "); below it",
    j, "fits, above it", length(s) - j, "\n")
cat("\nordstart's six flat points and the 5207 fit:\n")
six <- c("cs_sr_acat_seed4", "cs_sr_acat_seed1", "dpar_formula_seed4",
         "cum_sratio_seed14", "cum_sratio_seed18", "acat_cumcll_K7_seed6",
         "cum_sratio_seed16")
print(X[match(six, X$case), c("case", "ll", "code", "degen", "se_t", "se_b",
                              "maxpar")], row.names = FALSE)
cat("\nevery fit's largest coefficient SE above 20, with its threshold SE:\n")
print(X[(X$se_b > 20) %in% TRUE, c("case", "degen", "code", "se_t", "se_b",
                                   "maxpar")], row.names = FALSE)
