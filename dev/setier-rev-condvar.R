# Reviewer of lane setier: ranef(condVar = TRUE) against lme4 and
# glmmTMB, on a healthy fit and on fits that lost a standard error.
#   Rscript dev/setier-rev-condvar.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
has_tmb <- requireNamespace("glmmTMB", quietly = TRUE)
cat("arm", arm, find.package("frmtmb"), " glmmTMB", has_tmb, "\n")
q <- function(x) sprintf("%.4g", x)
cmp <- function(lab, ff, data, grp, family = gaussian()) {
  f <- suppressMessages(suppressWarnings(frm(ff, data = data,
                                              family = family)))
  rf <- suppressWarnings(ranef(f, condVar = TRUE))
  m <- suppressMessages(suppressWarnings(
    if (identical(family$family, "gaussian")) lmer(ff, data = data, REML = FALSE)
    else glmer(ff, data = data, family = family)))
  rl <- ranef(m, condVar = TRUE)
  cat("\n==", lab, "| lost:", paste(names(frmtmb:::sdr_of(f)$se_lost),
                                     collapse = ","), "\n")
  for (gname in grp) {
    a <- attr(rf[[gname]], "condSD")
    if (is.null(a)) a <- attr(rf[[gname]], "postVar")
    pl <- attr(rl[[gname]], "postVar")
    sl <- if (length(dim(pl)) == 3) sqrt(t(apply(pl, 3, diag))) else
      sqrt(pl)
    sl <- matrix(sl, nrow = nrow(rl[[gname]]))
    a <- matrix(a, nrow = nrow(rl[[gname]]))
    cat(" ", gname, ": frmtmb condSD / lme4 condSD, per column: median",
        paste(q(apply(a / sl, 2, median)), collapse = " "), " range",
        paste(q(range(a / sl)), collapse = " "), "\n")
    if (has_tmb) {
      gt <- suppressMessages(suppressWarnings(
        glmmTMB::glmmTMB(ff, data = data, family = family, REML = FALSE)))
      rt <- glmmTMB::ranef(gt, condVar = TRUE)$cond[[gname]]
      pt <- attr(rt, "condVar")
      st <- if (length(dim(pt)) == 3) sqrt(t(apply(pt, 3, diag))) else
        sqrt(pt)
      st <- matrix(st, nrow = nrow(rt))
      cat("   frmtmb / glmmTMB: median", paste(q(apply(a / st, 2, median)),
                                                collapse = " "),
          " range", paste(q(range(a / st)), collapse = " "), "\n")
    }
  }
}
s <- lme4::sleepstudy
cmp("healthy (Days | Subject)", Reaction ~ Days + (Days | Subject), s,
    "Subject")
s$a <- factor(s$Days %% 3)
cmp("boundary (1 | Subject/a)", Reaction ~ Days + (1 | Subject/a), s,
    c("Subject", "a:Subject"))
set.seed(30)
d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
cmp("boundary slope (1 + x | g), slope sd 0", y ~ x + (1 + x | g), d, "g")
cb <- lme4::cbpp
cmp("binomial cbpp (1 | herd)",
    cbind(incidence, size - incidence) ~ period + (1 | herd), cb, "herd",
    family = binomial())
