# Review of lane sampfront: the six watch callers give the same answers on
# ordinary draws, base against lane. Writes an RDS per arm; a second call
# with "compare" reads both and prints the largest difference per quantity.
#   Rscript dev/review-sampfront-same.R [base|lane|compare]
arm <- commandArgs(TRUE)[1]
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-sampfront/dev/review-sampfront-out"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(out_dir, "same-base.rds"))
  b <- readRDS(file.path(out_dir, "same-lane.rds"))
  for (k in union(names(a), names(b))) {
    x <- a[[k]]; y <- b[[k]]
    msg <- if (is.character(x) || is.character(y)) {
      paste("base:", substr(paste(x), 1, 90), "| lane:", substr(paste(y), 1, 90))
    } else if (!identical(dim(x), dim(y))) {
      paste("DIM base", paste(dim(x), collapse = "x"), "lane",
            paste(dim(y), collapse = "x"))
    } else {
      d <- max(abs(unlist(x) - unlist(y)), na.rm = TRUE)
      paste0("identical=", identical(x, y), " max|diff|=", format(d),
             " NA base/lane=", sum(is.na(unlist(x))), "/", sum(is.na(unlist(y))))
    }
    cat(sprintf("%-42s %s\n", k, msg))
  }
  quit(save = "no")
}
base <- identical(arm, "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, find.package("frmtmb.sample"), "\n")
lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
res <- list()
put <- function(k, expr) {
  set.seed(42)
  res[[k]] <<- tryCatch(suppressWarnings(suppressMessages(expr)),
                        error = function(e) paste("ERROR:", conditionMessage(e)))
}
set.seed(77)
G <- 8; n <- 12
dd <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n))
u <- rnorm(G, 0, 0.6)
dd$y <- 0.5 + 0.4 * dd$x + u[dd$g] + rnorm(G * n, 0, 0.7)
dd$cnt <- rpois(G * n, exp(0.2 + 0.3 * dd$x + u[dd$g]))
dd$ym <- ifelse(runif(G * n) < 0.4, rnorm(G * n, -2), rnorm(G * n, 2))
f1 <- q(bf(y ~ x + (1 | g)), data = dd)
p1 <- lap_pair(f1)
f2 <- q(bf(cnt ~ x + (1 | g)), family = poisson(), data = dd)
p2 <- lap_pair(f2)
fx <- q(bf(y ~ x + (0 + x | g)), data = dd)
px <- lap_pair(fx)
nd0 <- data.frame(x = c(0, 0, 0), g = factor(c(1, 2, 3), levels = 1:G))
for (w in c("full", "lap")) {
  put(paste0(w, " epred gaussian"), posterior_epred(p1[[w]]))
  put(paste0(w, " epred gaussian re=NA"),
      posterior_epred(p1[[w]], re_formula = NA))
  put(paste0(w, " linpred poisson"), posterior_linpred(p2[[w]]))
  put(paste0(w, " linpred poisson re=NA"),
      posterior_linpred(p2[[w]], re_formula = NA))
  put(paste0(w, " predict poisson"), posterior_predict(p2[[w]]))
  put(paste0(w, " predict poisson re=NA"),
      posterior_predict(p2[[w]], re_formula = NA))
  put(paste0(w, " hypothesis b_x"), hypothesis(p1[[w]], "x > 0")$samples)
  put(paste0(w, " hypothesis sd"),
      hypothesis(p1[[w]], "sd_g__Intercept > 0.1", class = NULL)$samples)
  put(paste0(w, " ce gaussian"),
      conditional_effects(p1[[w]], effects = "x")[[1]][, c("estimate__",
                                                         "lower__", "upper__")])
  put(paste0(w, " ce predict"),
      conditional_effects(p1[[w]], effects = "x",
                          method = "posterior_predict")[[1]][, c("estimate__",
                                                               "lower__")])
  put(paste0(w, " 1530 epred x=0"), posterior_epred(px[[w]], newdata = nd0))
  put(paste0(w, " 1530 predict x=0"), posterior_predict(px[[w]], newdata = nd0))
}
fm <- q(bf(ym ~ 1), family = mixture(gaussian(), gaussian()), data = dd,
        prior = c(set_prior("normal(-2, 1)", class = "Intercept", dpar = "mu1"),
                  set_prior("normal(2, 1)", class = "Intercept", dpar = "mu2")))
pm <- lap_pair(fm)
put("full pp_mixture", pp_mixture(pm$full, summary = FALSE))
fmv <- q(mvbf(bf(y ~ x), bf(cnt ~ x), rescor = FALSE),
         family = list(gaussian(), poisson()), data = dd)
pmv <- lap_pair(fmv)
put("full predict mv", posterior_predict(pmv$full))
fr <- q(mvbf(bf(y ~ x + (1 | g)), bf(ym ~ x), rescor = TRUE), data = dd)
pr <- lap_pair(fr)
put("full predict rescor", posterior_predict(pr$full))
put("lap predict rescor re=NA", posterior_predict(pr$lap, re_formula = NA))
saveRDS(res, file.path(out_dir, paste0("same-", arm, ".rds")))
cat("saved", length(res), "\n")
