# Claim 7: the ?frm claim that cs() "is not available ... in a
# multivariate fit". And claim 4's newdata-level question against brms,
# through validate_newdata, which needs no compile.
#   Rscript dev/csfactor-rev-docs.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
cat("== build", TAG, ": frmtmb", as.character(packageVersion("frmtmb")),
    " brms", as.character(packageVersion("brms")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 300L)
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  print(tryCatch(expr, error = function(e) paste("ERROR:", short(e))))
}

set.seed(77)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(d$fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n); u2 <- runif(n)
d$yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
d$yo2 <- ifelse(u2 < p1, 1L, ifelse(u2 < p1 + p2, 2L, 3L))

cat("\n########## does cs() work in a MULTIVARIATE fit? ##########\n")
mv <- tryCatch(suppressWarnings(suppressMessages(
        frm(frmtmb::bf(yo ~ x + cs(fc)) + frmtmb::bf(yo2 ~ z + cs(fc)),
            family = sratio(), data = d))),
      error = function(e) structure(list(m = short(e)), class = "revfail"))
if (inherits(mv, "revfail")) {
  cat("mv cs REFUSED: ", mv$m, "\n")
} else {
  cat("mv cs FITTED npar = ", length(mv$opt$par), " logLik = ",
      sprintf("%.7f", as.numeric(logLik(mv))), "\n", sep = "")
  go("fixef(mv)", rownames(fixef(mv)))
  go("variables(mv) bcs rows", grep("^bcs", variables(mv), value = TRUE))
  nd <- data.frame(x = 0, z = 0, fc = factor("c", levels = c("a","b","c")))
  go("fitted(mv, newdata at level c)", round(fitted(mv, newdata = nd), 4))
  # the univariate counterparts, to see whether the mv answer is the
  # per-response answer
  u1 <- frm(frmtmb::bf(yo ~ x + cs(fc)), family = sratio(), data = d)
  u2f <- frm(frmtmb::bf(yo2 ~ z + cs(fc)), family = sratio(), data = d)
  cat("logLik(mv) = ", sprintf("%.7f", as.numeric(logLik(mv))),
      "   sum of the two univariate = ",
      sprintf("%.7f", as.numeric(logLik(u1)) + as.numeric(logLik(u2f))),
      "\n", sep = "")
  cat("difference = ",
      sprintf("%.3e", as.numeric(logLik(mv)) -
                (as.numeric(logLik(u1)) + as.numeric(logLik(u2f)))),
      "\n", sep = "")
  go("frm_linpred(mv, resp = 'yo', newdata) at level c",
     round(frm_linpred(mv, newdata = nd, resp = "yo",
                       type = "response"), 4))
}

cat("\n########## a newdata level the fit never saw ##########\n")
ff <- frm(frmtmb::bf(yo ~ cs(fc)), family = sratio(), data = d)
ndbad <- data.frame(fc = factor("zz"))
go("frmtmb fitted() at an unseen level",
   fitted(ff, newdata = ndbad))
go("brms validate_newdata at an unseen level (no compile)",
   brms::standata(brms::bf(yo ~ cs(fc)), data = d,
                  family = brms::sratio(), newdata = ndbad))
go("brms validate_newdata at a level present in the fit",
   dim(brms::standata(brms::bf(yo ~ cs(fc)), data = d,
                      family = brms::sratio(),
                      newdata = data.frame(
                        fc = factor("c", levels = c("a","b","c")),
                        yo = 1L))$Xcs))

cat("\n########## variables() on a draws object: is the raw naming",
    "cs-specific? ##########\n")
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  library(frmtmb.sample)
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  ds1 <- tryCatch(frm_sample(frmtmb::bf(yo ~ x), family = sratio(),
                             data = d, chains = 1L, iter = 400L,
                             refresh = 0L, seed = 2L),
                  error = function(e) structure(list(m = short(e)),
                                                class = "revfail"))
  if (inherits(ds1, "revfail")) {
    cat("plain ordinal frm_sample FAILED: ", ds1$m, "\n")
  } else {
    go("variables(ds) for an ordinal fit with NO cs()", variables(ds1))
    go("fixef(ds) rows for the same", rownames(fixef(ds1)))
  }
}
cat("\nDONE ", TAG, "\n")
