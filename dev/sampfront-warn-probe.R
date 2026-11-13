# Lane sampfront fix round: which call of the cancel test warns.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
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
set.seed(77)
G <- 8; n <- 12
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u0 <- rnorm(G, 0, 0.8); u1 <- rnorm(G, 0, 0.8)
dd$y <- 0.5 + 0.4 * dd$x + u0[dd$g] + u1[dd$g] * dd$x + rnorm(nrow(dd), 0, 0.5)
pa <- lap_pair(suppressMessages(frm(bf(y ~ x + (1 + x | g)), data = dd)))
nda <- data.frame(x = -1, g = factor(c(1, 2), levels = 1:G))
g1 <- data.frame(g = factor(1, levels = 1:G))
w <- function(lab, expr) {
  ws <- character()
  r <- withCallingHandlers(tryCatch(expr, error = function(e) "refused"),
    warning = function(x) { ws <<- c(ws, conditionMessage(x))
      invokeRestart("muffleWarning") })
  cat(lab, ": ", if (identical(r, "refused")) "refused" else "ANSWERED",
      " warnings=", length(ws), " ", paste(unique(ws), collapse = "; "),
      "\n", sep = "")
}
w("epred", posterior_epred(pa$lap, newdata = nda))
w("linpred", posterior_linpred(pa$lap, newdata = nda))
w("predict", posterior_predict(pa$lap, newdata = nda))
w("ce", conditional_effects(pa$lap, effects = "x", re_formula = NULL,
                            int_conditions = list(x = -1), conditions = g1))
w("ce predict", conditional_effects(pa$lap, effects = "x",
                                    re_formula = NULL,
                                    method = "posterior_predict",
                                    int_conditions = list(x = -1),
                                    conditions = g1))
