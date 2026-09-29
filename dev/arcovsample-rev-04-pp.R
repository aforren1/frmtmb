# REVIEW script 04, claim 3: no behaviour change elsewhere.
#
#   Rscript dev/arcovsample-rev-04-pp.R <lane|ref> <out.rds>
#
# Both builds run the SAME code on the SAME deterministic draws objects
# (stanfit = NULL, a fixed parameter matrix, no sampler), so the only
# difference between the two runs is the installed library.
# posterior_predict() draws random numbers, so each call is wrapped in
# its own set.seed(); the fifth case is a LAPLACE draws object, made by
# dropping the random-effect columns, which is what draws_is_laplace()
# counts.

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]; out <- a[2L]
LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
newx <- all(c("arma_cond_resp", "arma_cond_dpars") %in%
              getNamespaceExports("frmtmb"))
cat("ARM ", arm, " frmtmb at ",
    dirname(system.file("DESCRIPTION", package = "frmtmb")),
    " newexports=", newx, "\n", sep = "")
stopifnot(identical(newx, identical(arm, "lane")))

set.seed(808L)
ng <- 5L
lens <- c(5L, 8L, 4L, 7L, 6L)
dd <- do.call(rbind, lapply(seq_len(ng), function(i) {
  data.frame(g = factor(i, levels = seq_len(ng)), t = seq_len(lens[i]))
}))
n <- nrow(dd)
dd$x <- rnorm(n)
dd$z <- runif(n, -1, 1)
dd$y <- 0.7 + 0.5 * dd$x + rnorm(n, 0, 0.8)
dd$y2 <- -0.2 + 0.3 * dd$x + rnorm(n, 0, 0.6)
dd <- dd[sample(n), ]
rownames(dd) <- NULL

# A draws matrix whose values are FIXED, so the two builds get the same
# numbers without a sampler. Every label but lp__ must be filled, or a
# zero sigma would quietly make the comparison trivial.
mk_draws <- function(fit, nd = 6L, laplace = FALSE) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(0, nd, length(lab), dimnames = list(NULL, lab))
  set.seed(4321L)
  for (j in seq_along(lab)) {
    nm <- lab[j]
    m[, j] <- if (nm == "lp__") {
      0
    } else if (grepl("^sigma", nm) || grepl("^sd_", nm)) {
      exp(rnorm(nd, log(0.8), 0.05))
    } else if (nm == "nu") {
      exp(rnorm(nd, log(8), 0.1))
    } else if (grepl("^thetaac", nm)) {
      rnorm(nd, 0.35, 0.05)
    } else if (grepl("^rescor", nm) || grepl("^cor_", nm)) {
      rnorm(nd, 0.3, 0.02)
    } else if (grepl("^thetar", nm)) {
      rnorm(nd, 0.3, 0.02)
    } else {
      rnorm(nd, 0.4, 0.1)
    }
  }
  if (laplace) {
    keep <- c(frmtmb::brms_par_labels(fit, include_random = FALSE), "lp__")
    m <- m[, intersect(colnames(m), keep), drop = FALSE]
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

fits <- list()
fits$condarma <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1)),
                     family = gaussian(), data = dd,
                     dry_run = "objective")
fits$covar <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
                  data = dd, dry_run = "objective")
fits$plain <- frm(bf(y ~ x), family = gaussian(), data = dd,
                  dry_run = "objective")
fits$rescor <- frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
                   family = gaussian(), data = dd, dry_run = "objective")
fits$re <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               dry_run = "objective")
fits$cosy <- frm(bf(y ~ x + cosy(t, g)), family = gaussian(), data = dd,
                 dry_run = "objective")
fits$unstr <- frm(bf(y ~ x + unstr(t, g)), family = gaussian(), data = dd,
                  dry_run = "objective")
dm <- dd; dm$y[c(3L, 12L)] <- NA
fits$mi <- frm(bf(y | mi() ~ x + ar(t, g)) + gaussian(), data = dm,
               dry_run = "objective")

res <- list(arm = arm, newexports = newx)

pp <- function(nm, laplace = FALSE) {
  ds <- mk_draws(fits[[nm]], laplace = laplace)
  set.seed(707L)
  p <- tryCatch(posterior_predict(ds), error = function(e)
    paste("ERR:", conditionMessage(e)))
  set.seed(707L)
  e <- tryCatch(posterior_epred(ds), error = function(e)
    paste("ERR:", conditionMessage(e)))
  list(predict = p, epred = e)
}
res$pp_condarma <- pp("condarma")
res$pp_covar <- pp("covar")
res$pp_plain <- pp("plain")
res$pp_rescor <- pp("rescor")
res$pp_laplace <- pp("re", laplace = TRUE)
res$pp_re <- pp("re")

msg <- function(nm, f) {
  ds <- mk_draws(fits[[nm]])
  conditionMessage(tryCatch({ f(ds); simpleError("NO ERROR") },
                            error = identity))
}
for (nm in c("covar", "cosy", "unstr", "mi", "condarma")) {
  res[[paste0("msg_loglik_", nm)]] <- msg(nm, function(d) log_lik(d))
  res[[paste0("msg_loo_", nm)]] <- msg(nm, function(d) loo(d))
}

saveRDS(res, out)
cat("wrote ", out, "\n", sep = "")
for (nm in grep("^msg_", names(res), value = TRUE)) {
  cat("\n[", nm, "]\n", res[[nm]], "\n", sep = "")
}
cat("\nshapes:\n")
for (nm in grep("^pp_", names(res), value = TRUE)) {
  v <- res[[nm]]$predict
  cat("  ", nm, ": ", if (is.character(v)) v else
    paste(dim(v), collapse = " x "), "\n", sep = "")
}
cat("DONE\n")
