# Reviewer probe 8: column order of a cross-response block in brms and
# frmtmb, which matters under prior = "flat" because flat on the
# Cholesky factor is not exchangeable. Also brms's program for lb on
# class cor. Seed 1.
arm <- commandArgs(trailingOnly = TRUE)[1L]
set.seed(1)
G <- 10L; n <- 6L
d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
d$y <- rnorm(G * n); d$y2 <- rnorm(G * n)
d$yn <- 2 * exp(0.3 * d$x) + rnorm(G * n, 0, 0.3)
fm <- list(
  mv = function(bf) bf(y ~ x + (1 + x | p | g)) + bf(y2 ~ x + (1 | p | g)),
  nl = function(bf) bf(yn ~ a * exp(k * x), a ~ 1 + (1 | q | g),
                       k ~ 1 + (1 | q | g), nl = TRUE))
if (identical(arm, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r14",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(brms))
  for (nm in names(fm)) {
    f <- fm[[nm]](brms::bf)
    if (nm == "mv") f <- f + set_rescor(FALSE)
    re <- brms:::frame_re(brms::brmsterms(f), d)
    cat(nm, ": ", paste(paste0(re$resp, re$nlpar, ":", re$coef),
                        collapse = " "), "\n", sep = "")
  }
  sc <- stancode(bf(y ~ x + (1 + x | g)), data = d,
                 prior = set_prior("lkj(2)", class = "cor", lb = 0))
  cat(grep("L_1|lkj", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  r <- tryCatch(rstan::stanc(model_code = sc), error = function(e) e)
  cat("stanc with lb = 0 on cor:", if (inherits(r, "error"))
    substr(conditionMessage(r), 1, 200) else "parses", "\n")
} else {
  source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
  for (nm in names(fm)) {
    f <- fm[[nm]](frmtmb::bf)
    if (nm == "mv") f <- f + set_rescor(FALSE)
    des <- frmtmb:::prior_design(f, d, NULL, list())
    for (bk in des$frame$re_blocks) {
      cat(nm, ": ", paste(bk$cnms, collapse = " "), "\n", sep = "")
    }
    fit <- suppressWarnings(suppressMessages(frm(f, data = d)))
    v <- grep("^cor_", rownames(summary(fit)$random[[1]] %||% matrix()), value = TRUE)
    nmz <- grep("^cor_", frmtmb:::brms_par_labels(fit, include_random = FALSE),
                value = TRUE)
    cat("  cor labels:", nmz, "\n")
  }
}
