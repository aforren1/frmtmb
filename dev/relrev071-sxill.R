# Release reviewer of 0.71.0: lane illse's factor route on lane
# famsimplex's families. A degree-6 raw polynomial on [1, 2] against the
# same model in poly()'s orthogonal basis: the fitted values and their
# standard errors on each category predictor must agree whatever the
# basis, so the raw fit's are checked against the orthogonal fit's.
#   Rscript dev/relrev071-sxill.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(5)
n <- 600
x <- runif(n, 1, 2)
e <- cbind(0, sin(3 * x), cos(2 * x) - 0.3)
P <- exp(e) / rowSums(exp(e))
G <- matrix(rgamma(3 * n, shape = 15 * P), n)
Y <- G / rowSums(G); colnames(Y) <- c("a", "b", "c")
yc <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
d <- data.frame(x = x, yc = yc)
d$Y <- Y
for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
raw <- stats::as.formula(paste("~", paste0("x", 1:6, collapse = " + ")))
nd <- d[seq(1, n, length.out = 25), ]
cmp <- function(lbl, fam, lhs) {
  fr <- suppressWarnings(frm(stats::as.formula(paste(lhs,
          paste0("x", 1:6, collapse = " + "), sep = " ~ ")),
          family = fam, data = d))
  fo <- suppressWarnings(frm(stats::as.formula(paste(lhs,
          "~ poly(x, 6)")), family = fam, data = d))
  cat(sprintf("%s: logLik raw %.6f orth %.6f | codes %d %d | route raw %s\n",
              lbl, as.numeric(logLik(fr)), as.numeric(logLik(fo)),
              fr$opt$convergence, fo$opt$convergence,
              !is.null(ns$sdr_of(fr)$cov_fixed_root)))
  l <- ns$sdr_of(fr)$se_lost
  cat("  raw lost:", if (length(l)) paste(names(l), l) else "-", "\n")
  for (dp in c("mub", "muc")) {
    pr <- tryCatch(frm_linpred(fr, newdata = nd, dpar = dp, se.fit = TRUE),
                   error = function(e) conditionMessage(e))
    po <- frm_linpred(fo, newdata = nd, dpar = dp, se.fit = TRUE)
    if (is.character(pr)) { cat("  ", dp, "ERROR", pr, "\n"); next }
    cat(sprintf("  %s: max |eta raw - orth| %.2e; max |se raw / se orth - 1| %.2e\n",
                dp, max(abs(pr$fit - po$fit)), max(abs(pr$se.fit / po$se.fit - 1))))
  }
}
cmp("dirichlet", dirichlet(), "Y")
cmp("logistic_normal", logistic_normal(), "Y")
cmp("categorical", categorical(), "yc")
