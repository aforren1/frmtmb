# Review of lane gpgap, item 4: non-centered at a RELATIVE nugget, in y's units.
# Absolute 1e-12 stopped at false convergence on 18 of 18 rows at y * 1000
# (ncp-60-A1000.txt). Also: does the centered 1e-6 reach its OWN closed-form
# maximum (seed 913 read short of it)?
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-ncp2.R <log> <n> <A>
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
a <- commandArgs(trailingOnly = TRUE)
n <- as.integer(a[1])
A <- as.numeric(a[2])
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
rel_corr <- function(nu) {
  function(theta, blk) {
    Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
    exp(-Q) + diag(nu, nrow(Q))
  }
}
abs_corr <- function(theta, blk) {
  Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
  exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
}
put <- function(val) {
  unlockBinding("gp_corr", ns); assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
closed <- function(x, y, s0, nug, st0 = NULL) {
  u <- sort(unique(x)); D2 <- outer(u, u, "-")^2; J <- match(x, u)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nug(s2)
    V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
  }
  best <- NULL
  for (st in c(list(c(0, log(A), 0, log(s0)), c(0, log(A) + 0.7, 0.3, log(s0)),
                    c(0, log(A) + 1.5, 0.5, log(s0) - 0.5)),
               if (!is.null(st0)) list(st0))) {
    o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
for (seed in c(911L, 913L)) {
  set.seed(seed)
  x <- round(runif(n, 0, 6), 3)
  e0 <- rnorm(n)
  for (s in c(1e-1, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)) {
    y <- A * (sin(2 * x) + s * e0)
    dd <- data.frame(x = x, y = y)
    fr <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = dd))
    bk <- Filter(function(k) identical(k$covstruct, "gp"), fr$frame$re_blocks)[[1]]
    th <- fr$estimates$theta[bk$theta_idx]
    own <- closed(x, y, A * s, function(s2) 1e-6 * s2)
    cl <- closed(x, y, A * s, function(s2) 1e-12 * max(A^2, 1) / max(A^2, 1))
    cat(sprintf(paste0("NC2 seed %d n %d A %g s %.0e | centered: ll %.4f sdgp %.4g conv %d; ",
                       "own closed form ll %.4f sdgp %.4g; centered short of it by %.4f | ",
                       "brms abs 1e-12 closed ll %.4f sdgp %.4g\n"),
                seed, n, A, s, -fr$opt$objective, exp(th[[1]]), fr$opt$convergence,
                -own$objective, exp(own$par[2]), fr$opt$objective - own$objective,
                -cl$objective, exp(cl$par[2])))
    for (v in list(list("rel 1e-8", rel_corr(1e-8)), list("rel 1e-10", rel_corr(1e-10)),
                   list("rel 1e-12", rel_corr(1e-12)), list("abs 1e-12", abs_corr))) {
      frn <- fr$frame
      frn$ncp_blocks <- 1L
      put(v[[2]])
      r <- tryCatch({
        fn <- frmtmb::build_objective(frn)
        tpl <- frn$par_template
        tpl$b[] <- 0
        t1 <- proc.time()[["elapsed"]]
        o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
        op <- suppressWarnings(nlminb(o$par, o$fn, o$gr,
                                      control = list(iter.max = 500, eval.max = 800)))
        h <- op$par[names(op$par) == "theta"]
        sprintf("conv %d gr %.1e ll %.4f sdgp %.4g %.1fs", op$convergence,
                max(abs(o$gr(op$par))), -op$objective, exp(h[[1]]),
                proc.time()[["elapsed"]] - t1)
      }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 60)))
      put(orig)
      cat(sprintf("   ncp %-10s %s\n", v[[1]], r))
    }
  }
}
