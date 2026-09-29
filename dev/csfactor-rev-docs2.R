# Follow-ups: is the mv fitted(newdata=) error cs-specific, what does
# brms do with an unseen newdata level, and is variables() on a draws
# object raw for the thresholds too or only for cs()?
#   Rscript dev/csfactor-rev-docs2.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
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
nd <- data.frame(x = 0, z = 0, fc = factor("c", levels = c("a","b","c")))

cat("\n### is the mv fitted(newdata=) failure cs-specific?\n")
mv_cs <- suppressWarnings(suppressMessages(
  frm(frmtmb::bf(yo ~ x + cs(fc)) + frmtmb::bf(yo2 ~ z + cs(fc)),
      family = sratio(), data = d)))
mv_no <- suppressWarnings(suppressMessages(
  frm(frmtmb::bf(yo ~ x + fc) + frmtmb::bf(yo2 ~ z + fc),
      family = sratio(), data = d)))
go("fitted(mv WITH cs, newdata)", dim(fitted(mv_cs, newdata = nd)))
go("fitted(mv WITHOUT cs, newdata)", dim(fitted(mv_no, newdata = nd)))
go("fitted(mv WITH cs, in sample)", dim(fitted(mv_cs)))
go("fitted(mv WITHOUT cs, in sample)", dim(fitted(mv_no)))
go("predict(mv WITH cs, newdata)",
   dim(predict(mv_cs, newdata = nd, ndraws = 50)))
go("posterior_epred(mv WITH cs, newdata)",
   dim(posterior_epred(mv_cs, newdata = nd, ndraws = 50)))

cat("\n### brms on an unseen newdata level, through validate_newdata\n")
bfit <- tryCatch(brms::brm(brms::bf(yo ~ cs(fc)), data = d,
                           family = brms::sratio(), empty = TRUE),
                 error = function(e) structure(list(m = short(e)),
                                               class = "revfail"))
if (inherits(bfit, "revfail")) {
  cat("brm(empty) failed: ", bfit$m, "\n")
} else {
  go("brms validate_newdata, level 'zz'",
     brms:::validate_newdata(data.frame(fc = factor("zz"), yo = 1L), bfit))
  go("brms validate_newdata, level 'c' (dim of the frame)",
     dim(brms:::validate_newdata(
       data.frame(fc = factor("c", levels = c("a","b","c")), yo = 1L),
       bfit)))
}

cat("\n### variables() on a draws object: cs-specific or not?\n")
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  library(frmtmb.sample)
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
  ds2 <- tryCatch(frm_sample(frmtmb::bf(yo ~ x), family = gaussian(),
                             data = d, chains = 1L, iter = 400L,
                             refresh = 0L, seed = 2L),
                  error = function(e) structure(list(m = short(e)),
                                                class = "revfail"))
  if (!inherits(ds2, "revfail")) {
    go("variables(ds) for a gaussian fit", variables(ds2))
  }
}
cat("\nDONE ", TAG, "\n")
