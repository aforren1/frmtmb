# Reviewer follow-ups to gpby-rev-core.R: gr = FALSE positions in row
# order, class "sd" on gp + (1 | g) read by name, the asymmetry of
# extra_cov at a shared 15-digit key, and timing done right (each round
# re-evaluates the call).
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d <- gpby_data()
tm <- function(cl, k = 5, env = parent.frame()) {
  min(vapply(seq_len(k), function(i) {
    system.time(eval(cl, env))[["elapsed"]]
  }, 0))
}
if (arm == "lane") {
  sdat <- suppressMessages(brms::standata(brms::bf(y ~ gp(x, by = f,
                                                          gr = FALSE)),
                                          data = d, family = gaussian()))
  fr <- frm(bf(y ~ gp(x, by = f, gr = FALSE)), data = d, dry_run = "frame")
  for (j in 1:3) {
    gi <- fr$linpreds[["y.mu"]]$gps[[j]]
    bk <- fr$re_blocks[[gi$block_id]]
    cat("gr = FALSE sub", j, "identical pos/div to Xgp:",
        identical(unname(gi$positions / bk$gp_lscale_div),
                  unname(as.matrix(sdat[[paste0("Xgp_1_", j)]]))), "\n")
  }
  d$g <- factor(rep(1:6, length.out = nrow(d)))
  vv <- function(f) {
    h <- frmtmb:::gp_brms_values(f, f$estimates$theta)
    s <- summary(f)$random
    c(h, sd_g = exp(f$estimates$theta[
      setdiff(seq_along(f$estimates$theta),
              unlist(lapply(Filter(function(b) b$covstruct == "gp",
                                   f$frame$re_blocks),
                            `[[`, "theta_idx")))]))
  }
  f0 <- frm(bf(y ~ gp(x) + (1 | g)), data = d)
  f1 <- frm(bf(y ~ gp(x) + (1 | g)), data = d,
            prior = set_prior("exponential(50)", class = "sd"))
  f2 <- frm(bf(y ~ gp(x) + (1 | g)), data = d,
            prior = set_prior("exponential(50)", class = "sdgp"))
  print(rbind(flat = vv(f0), class_sd = vv(f1), class_sdgp = vv(f2)))
  fit <- frm(bf(y ~ gp(x)), data = d)
  nd <- data.frame(x = c(7, 0.1 + 0.2 + 7 - 7, 0.3 + 7 - 7, 8))
  nd$x[2] <- 7.3; nd$x[3] <- 7.1 + 0.2
  cat("keys:", frmtmb:::pos_rowkey(as.matrix(nd$x)), "| x2 == x3:",
      nd$x[2] == nd$x[3], "\n")
  E <- frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)$extra_cov
  cat(sprintf("asymmetry max |E - t(E)| %.3e relative to max|E| %.3e\n",
              max(abs(E - t(E))), max(abs(E)) ))
}
fit <- frm(bf(y ~ gp(x)), data = d)
fs <- frm(bf(y ~ s(x, k = 10)), data = d)
big <- data.frame(x = seq(-1, 9, length.out = 2000))
bigs <- data.frame(x = seq(0.1, 5.9, length.out = 2000))
has_ec <- "extra_cov" %in% names(formals(frm_lp_basis))
cat(sprintf("TIME %s gp(x) 2000 rows frm_lp_basis %.3f s\n", arm,
            tm(quote(frm_lp_basis(fit, newdata = big)))))
if (has_ec) {
  cat(sprintf("TIME %s gp(x) 2000 rows frm_lp_basis(extra_cov) %.3f s\n",
              arm, tm(quote(frm_lp_basis(fit, newdata = big,
                                         extra_cov = TRUE)))))
}
cat(sprintf("TIME %s gp(x) 2000 rows frm_curve %.3f s\n", arm,
            tm(quote(frm_curve(fit, newdata = big, nsim = 1000, seed = 1)),
               3)))
cat(sprintf("TIME %s s(x) 2000 rows frm_curve %.3f s\n", arm,
            tm(quote(frm_curve(fs, newdata = bigs, nsim = 1000, seed = 1)),
               3)))
g0 <- gc(reset = TRUE)
cv <- frm_curve(fit, newdata = big, nsim = 1000, seed = 1)
g1 <- gc()
cat(sprintf("MEM %s gp(x) 2000 rows frm_curve max Vcells %.0f MB\n", arm,
            g1[2, 6]))
g0 <- gc(reset = TRUE)
cv <- frm_curve(fs, newdata = bigs, nsim = 1000, seed = 1)
g1 <- gc()
cat(sprintf("MEM %s s(x) 2000 rows frm_curve max Vcells %.0f MB\n", arm,
            g1[2, 6]))
cat("DONE\n")
