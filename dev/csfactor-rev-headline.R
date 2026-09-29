# Spot-check the headline numbers NEWS.md quotes: seed 405, n = 500,
# frm(bf(yo ~ cs(fc)), sratio()), logLik -446.9233 before and -446.6361
# after, and fitted() at factor("c") equal to level c's empirical shares.
#   Rscript dev/csfactor-rev-headline.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
# tests/testthat/test-cs-factor.R's csf_data(), verbatim
set.seed(405)
n <- 500
x <- stats::rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- stats::plogis(-0.3 + eff)
p2 <- (1 - p1) * stats::plogis(0.5 - eff)
u <- stats::runif(n)
d <- data.frame(x = x, fc = fc,
                yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fcb <- as.numeric(d$fc == "b"); d$fcc <- as.numeric(d$fc == "c")
ff <- frm(bf(yo ~ cs(fc)), family = sratio(), data = d)
fd <- frm(bf(yo ~ cs(fcb) + cs(fcc)), family = sratio(), data = d)
cat("logLik(cs(fc))       = ", sprintf("%.7f", as.numeric(logLik(ff))),
    "  npar = ", length(ff$opt$par), "\n", sep = "")
cat("logLik(hand dummies) = ", sprintf("%.7f", as.numeric(logLik(fd))),
    "  npar = ", length(fd$opt$par), "\n", sep = "")
cat("fixef rows: ", paste(rownames(fixef(ff)), collapse = ", "), "\n",
    sep = "")
P <- frm_linpred(ff, newdata = data.frame(fc = factor("c")),
                 type = "response")
cat("fitted at factor('c') = ", paste(sprintf("%.4f", P), collapse = ", "),
    "\n", sep = "")
E <- prop.table(table(d$fc, d$yo), 1L)
cat("empirical shares, level c = ",
    paste(sprintf("%.4f", E["c", ]), collapse = ", "), "\n", sep = "")
cat("max abs difference from the empirical shares = ",
    sprintf("%.3e", max(abs(as.numeric(P) - as.numeric(E["c", ])))),
    "\n", sep = "")
Pa <- frm_linpred(ff, newdata = data.frame(
  fc = factor("a", levels = c("a", "b", "c"))), type = "response")
cat("fitted at level a = ", paste(sprintf("%.4f", Pa), collapse = ", "),
    "\n", sep = "")
# in sample, against the hand-dummy fit: bitwise or not
F1 <- frm_linpred(ff, type = "response")
F2 <- frm_linpred(fd, type = "response")
cat("in-sample max abs difference from the hand-dummy fit = ",
    sprintf("%.17g", max(abs(F1 - F2))), "  identical = ",
    identical(F1, F2), "\n", sep = "")
cat("\nDONE ", TAG, "\n")
