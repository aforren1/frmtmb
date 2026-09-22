# Reviewer, priority 1: an OUTSIDE judge for the standard errors the
# interop seams report. frmtmb (base and lane) against the same model
# in glmmTMB, in lme4, in brms (posterior sd of the same quantity) and
# against a parametric bootstrap computed here without any frmtmb
# covariance on its path.
#
#   Rscript dev/shapes-rev-judge.R lane|base

arg <- commandArgs(trailingOnly = TRUE)
which_lib <- if (identical(arg[1], "base")) "base" else "lane"
lib <- if (which_lib == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/shapes-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
B <- readRDS(file.path(TREE, "dev/shapes-rev-brmsref.rds"))
dd <- rev_data()
cat("lib:", lib, "\n")

sl <- function(m) {
  s <- as.data.frame(suppressWarnings(marginaleffects::avg_slopes(m)))
  data.frame(term = as.character(s$term),
             contrast = if ("contrast" %in% names(s))
               as.character(s$contrast) else "",
             est = s$estimate, se = s$std.error,
             stringsAsFactors = FALSE)
}
show <- function(lab, x) {
  cat("  ", lab, "\n", sep = "")
  if (inherits(x, "try-error") || is.null(x)) { cat("    n/a\n"); return() }
  for (i in seq_len(nrow(x)))
    cat(sprintf("    %-10s %-10s est %10.5f  se %10.5f\n", x$term[i],
                x$contrast[i], x$est[i], x$se[i]))
}

# marginaleffects cannot run on a brmsfit here (it needs `collapse`,
# which is not installed and the user library is read only), so the
# brms judge is its POSTERIOR SD of the coefficient itself. On an
# identity link the average marginal effect of x IS b_x and a factor
# level's contrast IS its coefficient, so the two questions coincide
# exactly; on a logit link they do not, and the parametric bootstrap
# at the end is the judge there instead.
brms_fx <- function(nm) {
  fx <- B[[nm]]$fixef
  if (is.null(fx) || inherits(fx, "revErr")) { cat("   brms n/a\n"); return() }
  cat("   brms fixef (Estimate, posterior sd)\n")
  for (i in seq_len(nrow(fx)))
    cat(sprintf("    %-16s est %10.5f  se %10.5f\n", rownames(fx)[i],
                fx[i, "Estimate"], fx[i, "Est.Error"]))
}

cat("\n############ gaussian: y ~ x + f ############\n")
fg <- frm(bf(y ~ x + f) + gaussian(), data = dd)
show("frmtmb", try(sl(fg), silent = TRUE))
show("lm", try(sl(stats::lm(y ~ x + f, data = dd)), silent = TRUE))
show("glmmTMB", try(sl(glmmTMB::glmmTMB(y ~ x + f, data = dd)),
                    silent = TRUE))
brms_fx("gaussian")

cat("\n############ binomial: bin ~ x + z ############\n")
fb <- frm(bf(bin ~ x + z) + bernoulli(), data = dd)
show("frmtmb", try(sl(fb), silent = TRUE))
show("glm", try(sl(stats::glm(bin ~ x + z, binomial(), data = dd)),
                silent = TRUE))
show("glmmTMB", try(sl(glmmTMB::glmmTMB(bin ~ x + z, family = binomial(),
                                        data = dd)), silent = TRUE))
brms_fx("binomial")

cat("\n############ mixed: ymix ~ x + (1|g) ############\n")
fm <- frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd)
show("frmtmb", try(sl(fm), silent = TRUE))
show("lme4", try(sl(lme4::lmer(ymix ~ x + (1 | g), data = dd)),
                 silent = TRUE))
show("glmmTMB", try(sl(glmmTMB::glmmTMB(ymix ~ x + (1 | g), data = dd)),
                    silent = TRUE))
brms_fx("mixed")

cat("\n############ ordinal: ord ~ x, cumulative ############\n")
fo <- frm(bf(ord ~ x) + cumulative(), data = dd)
show("frmtmb", try(sl(fo), silent = TRUE))
show("MASS::polr", try(sl(MASS::polr(ord ~ x, data = dd, Hess = TRUE)),
                       silent = TRUE))
brms_fx("ordinal")

cat("\n############ an INDEPENDENT parametric bootstrap ############\n")
cat("# the binomial average marginal effect of x, resampled from the\n")
cat("# fitted glm and refitted, with no frmtmb covariance on its path\n")
set.seed(31415)
g0 <- stats::glm(bin ~ x + z, binomial(), data = dd)
Rb <- 400L
amex <- numeric(Rb); amez <- numeric(Rb)
for (b in seq_len(Rb)) {
  db <- dd
  db$bin <- stats::rbinom(nrow(dd), 1, stats::fitted(g0))
  gb <- suppressWarnings(stats::glm(bin ~ x + z, binomial(), data = db))
  cf <- stats::coef(gb)
  eta <- cf[1] + cf[2] * dd$x + cf[3] * dd$z
  w <- stats::dlogis(eta)
  amex[b] <- mean(w * cf[2]); amez[b] <- mean(w * cf[3])
}
cat(sprintf("    bootstrap R = %d:  x se %10.5f   z se %10.5f\n", Rb,
            stats::sd(amex), stats::sd(amez)))
