# Lane gpfield: shared helpers. Arm "lane" reads wt-gpfield-lib over
# rellib-r21; arm "base" reads rellib-r21 alone (frmtmb 0.82.0 at ec8178d3).
GF_WT <- "C:/Users/adf44/source/r/frmtmb-wt-rel083"
GF_LIB <- "C:/Users/adf44/source/r/rellib-r22"
GF_BASE <- "C:/Users/adf44/source/r/rellib-r21"
GF_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
GF_OUT <- file.path(GF_WT, "dev/rel083-log/kt/gpfield")
dir.create(GF_OUT, showWarnings = FALSE)
gf_libs <- function(arm = c("lane", "base")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") GF_LIB else GF_BASE, GF_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(GF_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "arm", arm, "\n")
  invisible(arm)
}
# A fit with its warnings kept (class and text), never an error thrown.
gf_fit <- function(...) {
  w <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(e) {
    w <<- c(w, paste0(class(e)[1], ": ", conditionMessage(e)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  list(fit = f, warn = w)
}
gf_put <- function(name, f) {
  ns <- asNamespace("frmtmb")
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
# median nearest-neighbour gap of distinct 1-D positions
gf_mednn <- function(x) {
  u <- sort(unique(x))
  g <- diff(u)
  stats::median(pmin(c(Inf, g), c(g, Inf)))
}
# smallest eigenvalue of the inner Hessian at outer point par (dense)
gf_inner_ev <- function(obj, par) {
  env <- obj$env
  if (!length(env$random)) return(NA_real_)
  obj$fn(par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  min(eigen(0.5 * (H + t(H)), symmetric = TRUE, only.values = TRUE)$values)
}
# marginal by importance sampling (lane gpstop's gpstop-stis2.R is_ll():
# multivariate t, df 4, proposal precision the inner Hessian's eigenvalues
# floored at 1, set.seed(9001))
gf_is <- function(obj, par, ndraw = 20000L, df = 4, seed = 9001L) {
  env <- obj$env
  lap <- -obj$fn(par)
  full <- env$last.par
  r <- env$random
  H <- as.matrix(env$spHess(full, random = TRUE))
  H <- 0.5 * (H + t(H))
  e <- eigen(H, symmetric = TRUE)
  ev <- e$values
  P <- e$vectors %*% diag(pmax(ev, 1)) %*% t(e$vectors)
  u0 <- full[r]
  R <- chol(0.5 * (P + t(P)))
  k <- length(r)
  set.seed(seed)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
  dev <- t(backsolve(R, t(z))) / s
  lq <- lgamma((df + k) / 2) - lgamma(df / 2) - k / 2 * log(df * pi) +
    sum(log(diag(R))) - (df + k) / 2 * log1p(rowSums(z^2) / s^2 / df)
  lj <- vapply(seq_len(ndraw), function(i) {
    x <- full
    x[r] <- u0 + dev[i, ]
    -env$f(x, order = 0)
  }, 0)
  lw <- lj - lq
  m <- max(lw)
  w <- exp(lw - m)
  h <- seq_len(ndraw) <= ndraw / 2
  c(laplace = lap, is = m + log(mean(w)), h1 = m + log(mean(w[h])),
    h2 = m + log(mean(w[!h])), ess = sum(w)^2 / sum(w^2),
    ev_min = min(ev))
}
# the named functions of a source file of this worktree, put into the
# loaded namespace (a check of an edit before the library can be
# reinstalled: a sweep is reading it)
gf_patch <- function(file, names) {
  ex <- parse(file.path(GF_WT, file), keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% names) {
      f <- eval(e[[3]], asNamespace("frmtmb"))
      gf_put(as.character(e[[2]]), f)
      cat("PATCH", as.character(e[[2]]), "\n")
    }
  }
}
