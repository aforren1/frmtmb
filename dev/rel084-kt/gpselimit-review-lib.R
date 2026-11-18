# rel084-kt copy of dev/reviews/gpselimit-review-lib.R: RV_WT = the rel084 tree (stan cache); RV_OUT = dev/rel084-log/kt/gpselimit-out; arm lane = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Review of lane gpselimit: shared helpers. Arm "lane" reads
# wt-gpselimit-lib over rellib-r22; arm "base" reads rellib-r22 alone.
RV_WT <- "C:/Users/adf44/source/r/frmtmb-wt-rel084"
RV_OUT <- file.path(RV_WT, "dev/rel084-log/kt/gpselimit-out")
dir.create(RV_OUT, showWarnings = FALSE)
rv_libs <- function(arm) {
  base <- "C:/Users/adf44/source/r/rellib-r22"
  user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
  lp <- switch(arm,
               lane = "C:/Users/adf44/source/r/rellib-r23",
               base = base, stop("arm"))
  .libPaths(c(lp, user))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(RV_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "arm", arm, "\n")
}
# A fit with every condition kept (all classes and text); errors returned.
rv_fit <- function(expr) {
  conds <- list()
  f <- tryCatch(withCallingHandlers(expr, warning = function(e) {
    conds[[length(conds) + 1L]] <<- list(type = "warning",
                                         class = class(e),
                                         text = conditionMessage(e))
    invokeRestart("muffleWarning")
  }, message = function(e) {
    conds[[length(conds) + 1L]] <<- list(type = "message",
                                         class = class(e),
                                         text = conditionMessage(e))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  list(fit = f, conds = conds)
}
rv_cstr <- function(conds) {
  if (!length(conds)) return("none")
  paste(vapply(conds, function(c) {
    paste0(substr(c$type, 1, 1), ":", c$class[1],
           if (length(c$class) > 3) paste0("+", c$class[2]) else "")
  }, ""), collapse = ";")
}
# Quietly evaluate, keeping the conditions' first class.
rv_q <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, paste0("w:", class(e)[1]))
    invokeRestart("muffleWarning")
  }, message = function(e) {
    w <<- c(w, paste0("m:", class(e)[1]))
    invokeRestart("muffleMessage")
  })
  attr(v, "conds") <- w
  v
}
rv_sweep_row <- function(row, jitter = 0) {
  ks <- c(5, 8, 10, 15, 20, 30)
  sets <- list(c(n = 60, ls0 = 0.7), c(n = 100, ls0 = 1.5),
               c(n = 100, ls0 = 3), c(n = 200, ls0 = 1),
               c(n = 50, ls0 = 2))
  si <- (row - 1L) %/% length(ks) + 1L
  k <- ks[(row - 1L) %% length(ks) + 1L]
  n <- sets[[si]][["n"]]; ls0 <- sets[[si]][["ls0"]]
  set.seed(9600 + row)
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * ls0^2)) + diag(1e-8, n)
  f <- as.vector(t(chol(K)) %*% rnorm(n))
  y <- f + 0.3 * rnorm(n)
  if (jitter) y <- y * (1 + 2e-15 * jitter)
  list(d = data.frame(x = x, y = y), k = k, n = n, ls0 = ls0)
}
# The objective along a white block's ridges at moves t (positive toward
# lscale 0), minus the fit's, through the package's own probe function.
rv_dirs <- function(bk) {
  nl <- length(bk[["theta_idx"]]) - 1L
  if (nl == 1L) return(list(c((bk[["gp_D"]] %||% 1) / 2, -1)))
  lapply(seq_len(nl), function(j) {
    v <- numeric(nl + 1L); v[c(1L, 1L + j)] <- c(0.5, -1); v
  })
}
rv_ridge2 <- function(fit, bk, ts) {
  ns <- asNamespace("frmtmb")
  om <- ns$outer_par_map(fit)
  thpos <- which(om$comp == "theta")[bk$theta_idx]
  p <- fit$opt$par
  saved <- ns$obj_state_save(fit$obj)
  on.exit(ns$obj_state_restore(fit$obj, saved))
  pf <- ns$se_probe_fn(fit)
  f0 <- pf(p)
  lapply(rv_dirs(bk), function(v) vapply(ts, function(t) {
    q <- p; q[thpos] <- q[thpos] + t * v
    tryCatch(pf(q) - f0, error = function(e) NA_real_)
  }, 0))
}
rv_hsgp <- function(fit) {
  Filter(function(b) identical(b$covstruct, "hsgp"), fit$frame$re_blocks)
}
rv_ses <- function(fit) {
  s <- rv_q(frmtmb:::sdr_of(fit))
  se <- sqrt(pmax(diag(s$cov.fixed), 0))
  se[is.na(diag(s$cov.fixed))] <- NaN
  list(se = se, lost = s$se_lost, conds = attr(s, "conds"))
}
