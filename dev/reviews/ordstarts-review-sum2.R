# Review of lane ordstarts: the corrected count of section 3 ("silent
# short of a clean maximum"), read from the lane's logs plus the
# review's summary() run (dev/reviews/ordstarts-review-swap-*.tsv).
# Clean B: code 0, NO warning of any kind at the fit (the lane asked only
# for no degenerate warning), threshold SE <= 100, largest finite SE <=
# 100, holds under the lane's jitter. Silent A: no warning at the fit and
# none at summary() (the standard-error check of A often waits).
#   Rscript dev/reviews/ordstarts-review-sum2.R
O <- "dev/ordstarts-out/"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "")
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
X$key <- paste(X$family, X$re, X$n, X$seed)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ][match(A$key, X$key[X$arm == "starts10"]), ]
d <- B$ll - A$ll
up <- d > 1e-3 & !is.na(d)
R <- rbind(rd(paste0(O, "real.tsv")), rd(paste0(O, "pilot-refits.tsv")))
R$key <- paste(R$family, R$re, R$n, R$seed)
holds <- vapply(A$key, function(k) {
  r <- R[R$key == k, ]
  if (!nrow(r)) return(NA)
  f0 <- r$ll[r$check == "fit" & r$variant == "B"][1]
  jj <- r[r$check == "real" & r$variant == "B" & r$code %in% 0, ]
  !nrow(jj) || max(jj$ll) - f0 <= 1e-3
}, NA)
top1 <- vapply(A$key, function(k) {
  s <- R$note[R$key == k & R$check == "se" & R$variant == "B"][1]
  if (is.na(s)) return(NA_real_)
  suppressWarnings(as.numeric(sub(" at .*", "", sub("^[^ ]+ se ", "",
                                                     sub(";.*", "", s)))))
}, 0)
tseB <- suppressWarnings(pmax(B$tse1, B$tse2, na.rm = TRUE))
S <- do.call(rbind, lapply(Sys.glob("dev/reviews/ordstarts-review-swap-c*.tsv"),
                           rd))
S$key <- paste(S$family, "int", S$n, S$seed)
summw <- S$summwarnA[match(A$key, S$key)]
lane24 <- up & A$nwarn == 0 & B$code %in% 0 & !B$degen_warn &
  !(tseB > 100) & holds %in% TRUE & !(top1 > 100)
cat("lane's 24 reproduced:", sum(lane24), "\n")
cat("  of those, B warns at the fit:", sum(lane24 & B$nwarn > 0), ":",
    paste(A$key[lane24 & B$nwarn > 0], collapse = "; "), "\n")
cat("  of those, A warns at summary():", sum(lane24 & summw %in% 1:99),
    ":", paste(A$key[lane24 & summw %in% 1:99], collapse = "; "), "\n")
strict <- lane24 & B$nwarn == 0 & !(summw %in% 1:99)
cat("strict count (A silent through summary(), B no warning):",
    sum(strict), "; gain", signif(range(d[strict]), 3), "\n")
print(table(sub(" [0-9]+$", "", A$key[strict])))
cat("  of those, B at a group sd of zero (singular message):",
    sum(strict & B$nlost > 0 & grepl("Boundary \\(singular\\)", B$msgs)),
    "\n")
cat("\nA silent at fit but warns at summary(), (1 | g):",
    sum(A$nwarn == 0 & summw %in% 1:99, na.rm = TRUE), "of",
    sum(A$re == "int" & !is.na(summw)), "; on up moves",
    sum(up & A$nwarn == 0 & summw %in% 1:99), ":",
    paste(A$key[up & A$nwarn == 0 & summw %in% 1:99], collapse = "; "),
    "\n")
cat("\nhint fits (A degenerate):", sum(A$degen_warn %in% TRUE),
    "; B without the degenerate warning:",
    sum(A$degen_warn %in% TRUE & !(B$degen_warn %in% TRUE)), "\n")
cat("B gains the degenerate warning at the SAME logLik (|d| <= 1e-6):",
    sum(!(A$degen_warn %in% TRUE) & B$degen_warn %in% TRUE &
          abs(d) <= 1e-6, na.rm = TRUE), ":",
    paste(A$key[!(A$degen_warn %in% TRUE) & B$degen_warn %in% TRUE &
                  abs(d) <= 1e-6 & !is.na(d)], collapse = "; "), "\n")
