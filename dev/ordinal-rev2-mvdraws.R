# Reviewer re-check (punch round 1), lane ordinal: frmtmb.sample draws of
# multivariate models with ordinal responses. For every ordinal response,
# posterior_epred(resp =) at three draws against brms's R-side density at
# that draw's STORED columns, and every delta column against the spacing
# of its stored thresholds. Data seed 20261012, sampler seed 5.
# Output: dev/ordinal-rev2-log-mvdraws.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", find.package("frmtmb"), " frmtmb.sample",
    find.package("frmtmb.sample"), "\n")
set.seed(20261012)
n <- 240
d <- data.frame(x = rnorm(n), h = factor(sample(c("p", "q"), n, TRUE)))
cuts <- function(k) {
  u <- stats::rlogis(n) + 0.8 * d$x
  1L + rowSums(outer(u, seq(-1.5, 1.5, length.out = k - 1), ">"))
}
d$y <- cuts(5); d$y2 <- cuts(3); d$y3 <- cuts(4)
d$yg <- 1 + d$x + rnorm(n)

brms_P <- function(fam, link, eta, th) {
  get(paste0("d", fam), asNamespace("brms"))(seq_len(length(th) + 1L),
    eta = eta, thres = matrix(th, length(eta), length(th), byrow = TRUE),
    disc = 1, link = link)
}
check <- function(lab, f, fams, ords) {
  cat("\n==", lab, "==\n")
  tryCatch({
    fit <- frm(f, data = d, family = fams)
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
    v <- variables(ds)
    cat("  variables:", v[!grepl("^(lp__|lprior)", v)], "\n")
    M <- as.matrix(ds)
    for (o in ords) {
      r <- o$resp
      ep <- posterior_epred(ds, resp = r)
      worst <- 0
      for (i in c(1L, 77L, nrow(M))) {
        eta <- d$x * M[i, paste0("b_", r, "_x")]
        tc <- grep(paste0("^b_", r, "_Intercept[[]"), colnames(M), value = TRUE)
        if (isTRUE(o$gr)) {
          P <- matrix(NA_real_, n, dim(ep)[3])
          for (g in levels(d$h)) {
            cols <- grep(paste0("[[]", g, ","), tc, value = TRUE)
            rows <- which(d$h == g)
            P[rows, seq_len(length(cols) + 1L)] <-
              brms_P(o$fam, o$link, eta[rows], M[i, cols])
          }
        } else {
          P <- brms_P(o$fam, o$link, eta, M[i, tc])
        }
        worst <- max(worst, max(abs(ep[i, , ] - P), na.rm = TRUE))
      }
      cat(sprintf("  resp %s: max |epred - brms at stored columns| %.3g\n",
                  r, worst))
      for (dn in grep(paste0("^delta_", r, "(_[0-9]+)?$"), colnames(M),
                      value = TRUE)) {
        k <- sub(paste0("^delta_", r, "_?"), "", dn)
        tcs <- grep(paste0("^b_", r, "_Intercept[[]"), colnames(M), value = TRUE)
        if (nzchar(k)) tcs <- grep(paste0("[[]", levels(d$h)[as.integer(k)], ","),
                                   tcs, value = TRUE)
        sp <- M[, tcs[2]] - M[, tcs[1]]
        cat(sprintf("  %s: max |delta - spacing| %.3g (mean delta %.4f)\n",
                    dn, max(abs(M[, dn] - sp)), mean(M[, dn])))
      }
      if (isTRUE(o$stz)) {
        tcs <- grep(paste0("^b_", r, "_Intercept[[]"), colnames(M), value = TRUE)
        cat(sprintf("  resp %s: max |row sum of thresholds| %.3g\n", r,
                    max(abs(rowSums(M[, tcs, drop = FALSE])))))
      }
    }
    if ("yg" %in% names(fit$spec$responses)) {
      ep <- posterior_epred(ds, resp = "yg")
      ref <- M[, "b_yg_Intercept"] + outer(M[, "b_yg_x"], d$x)
      cat(sprintf("  resp yg: max |epred - (a + b x)| %.3g\n",
                  max(abs(ep - ref))))
    }
    ll <- log_lik(ds)
    cat("  log_lik dim:", dim(ll), " all finite:", all(is.finite(ll)), "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}

check("three ordinal responses: cumulative equidistant, sratio, acat probit sum_to_zero",
      bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x),
      list(cumulative(threshold = "equidistant"), sratio(),
           acat("probit", threshold = "sum_to_zero")),
      list(list(resp = "y", fam = "cumulative", link = "logit"),
           list(resp = "y2", fam = "sratio", link = "logit"),
           list(resp = "y3", fam = "acat", link = "probit", stz = TRUE)))
check("three ordinal responses, equidistant on all (ordered and not)",
      bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x),
      list(sratio(threshold = "equidistant"),
           cumulative(threshold = "equidistant"),
           cratio(threshold = "equidistant")),
      list(list(resp = "y", fam = "sratio", link = "logit"),
           list(resp = "y2", fam = "cumulative", link = "logit"),
           list(resp = "y3", fam = "cratio", link = "logit")))
check("ordinal plus gaussian: cumulative equidistant + gaussian",
      bf(y ~ x) + bf(yg ~ x), list(cumulative(threshold = "equidistant"),
                                   gaussian()),
      list(list(resp = "y", fam = "cumulative", link = "logit")))
check("gaussian first, then two ordinal (sum_to_zero, equidistant)",
      bf(yg ~ x) + bf(y3 ~ x) + bf(y ~ x),
      list(gaussian(), cumulative(threshold = "sum_to_zero"),
           acat(threshold = "equidistant")),
      list(list(resp = "y3", fam = "cumulative", link = "logit", stz = TRUE),
           list(resp = "y", fam = "acat", link = "logit")))
check("thres(gr = h) inside mv: cumulative equidistant grouped + sratio sum_to_zero",
      bf(y | thres(gr = h) ~ x) + bf(y2 ~ x),
      list(cumulative(threshold = "equidistant"),
           sratio(threshold = "sum_to_zero")),
      list(list(resp = "y", fam = "cumulative", link = "logit", gr = TRUE),
           list(resp = "y2", fam = "sratio", link = "logit", stz = TRUE)))
check("thres(gr = h) inside mv: flexible grouped sratio + cumulative",
      bf(y | thres(gr = h) ~ x) + bf(y3 ~ x),
      list(sratio(), cumulative()),
      list(list(resp = "y", fam = "sratio", link = "logit", gr = TRUE),
           list(resp = "y3", fam = "cumulative", link = "logit")))
