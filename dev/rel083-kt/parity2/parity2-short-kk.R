# Short check of parity2's fix round, F3: a numeric variable that is
# both a fixed predictor and a grouping factor (y ~ kk * x + (1 | kk)),
# at re_formula = NA (both packages' default) and NULL. Arm "brms":
# brms 2.23.0's own steps (update_re_terms(), get_all_effects(),
# prepare_cond_data()) on brm(empty = TRUE); arm "frmtmb": the lane.
#   Rscript dev/reviews/parity2-short-kk.R brms|frmtmb > <out> 2>&1
arm <- commandArgs(TRUE)[1]
set.seed(21)
n <- 200
e <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                h = factor(sample(c("u", "v"), n, TRUE)),
                grp = factor(sample(letters[1:8], n, TRUE)),
                kk = sample(1:6, n, TRUE))
e$y <- 1 + e$x + rnorm(6)[e$kk] + rnorm(n)
fo <- y ~ kk * x + (1 | kk)
if (identical(arm, "brms")) {
  .libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
  suppressMessages(library(brms))
  fit <- brm(fo, data = e, empty = TRUE)
  for (ref in list(NA, NULL)) {
    bt <- brms:::brmsterms(brms:::update_re_terms(fit$formula,
                                                  re_formula = ref))
    rsv <- brms:::rsv_vars(bt)
    gv <- brms:::get_group_vars(bt)
    nm <- vapply(brms:::get_all_effects(bt, rsv_vars = rsv), function(ev) {
      cd <- brms:::prepare_cond_data(
        model.frame(fit)[, ev, drop = FALSE],
        conditions = data.frame(row.names = 1),
        int_vars = brms:::get_int_vars(bt), group_vars = gv,
        resolution = 5, reorder = TRUE)
      paste(attr(cd, "effects"), collapse = ":")
    }, "")
    cat("brms re_formula =", deparse(ref), "group_vars:",
        paste(gv, collapse = " "), "names:", nm, "\n")
  }
} else {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r22",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", find.package("frmtmb"), "\n")
  fit <- suppressMessages(frm(bf(fo), family = gaussian(), data = e))
  for (ref in list(NA, NULL)) {
    ce <- suppressMessages(conditional_effects(fit, resolution = 5,
                                               re_formula = ref))
    cat("frmtmb re_formula =", deparse(ref), "names:", names(ce), "\n")
  }
}
cat("done\n")
