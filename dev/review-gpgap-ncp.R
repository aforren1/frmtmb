# Review of lane gpgap, item 4: the non-centered option at brms's nugget, on
# seeds that are not the lane's, against the closed form at that nugget and
# against the centered 1e-6. Also: can the kernel be factored at all at
# brms's absolute 1e-12 (Stan rejects a proposal, the tape returns NaN), and
# how does that scale with y's units?
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-ncp.R <log> <n> <A>
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
a <- commandArgs(trailingOnly = TRUE)
n <- as.integer(a[1])
A <- as.numeric(a[2])
ns <- asNamespace("frmtmb")
orig <- get("gp_corr", envir = ns)
abs_corr <- function(theta, blk) {
  Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
  exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
}
put <- function(val) {
  unlockBinding("gp_corr", ns); assign("gp_corr", val, envir = ns)
  lockBinding("gp_corr", ns)
}
closed <- function(x, y, s0) {
  u <- sort(unique(x)); D2 <- outer(u, u, "-")^2; J <- match(x, u)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + 1e-12
    V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
  }
  best <- NULL
  for (st in list(c(0, log(A), 0, log(s0)), c(0, log(A) + 0.7, 0.3, log(s0)),
                  c(0, log(A) + 1.5, 0.5, log(s0) - 0.5))) {
    o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
for (seed in c(911L, 912L, 913L)) {
  set.seed(seed)
  x <- round(runif(n, 0, 6), 3)
  e0 <- rnorm(n)
  for (s in c(1e-1, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)) {
    y <- A * (sin(2 * x) + s * e0)
    dd <- data.frame(x = x, y = y)
    t0 <- proc.time()[["elapsed"]]
    fr <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = dd))
    t_c <- proc.time()[["elapsed"]] - t0
    cl <- closed(x, y, A * s)
    ncpf <- function() {
      frn <- fr$frame
      frn$ncp_blocks <- 1L
      put(abs_corr)
      on.exit(put(orig))
      fn <- frmtmb::build_objective(frn)
      tpl <- frn$par_template
      tpl$b[] <- 0
      t1 <- proc.time()[["elapsed"]]
      o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
      op <- suppressWarnings(nlminb(o$par, o$fn, o$gr,
                                    control = list(iter.max = 500, eval.max = 800)))
      th <- op$par[names(op$par) == "theta"]
      sprintf("conv %d (%s) gr %.1e ll %.4f sdgp %.4g %.1fs", op$convergence,
              substr(op$message, 1, 16), max(abs(o$gr(op$par))), -op$objective,
              exp(th[[1]]), proc.time()[["elapsed"]] - t1)
    }
    r <- tryCatch(ncpf(), error = function(e) {
      paste("ERROR", substr(conditionMessage(e), 1, 70))
    })
    sdgp_cl <- exp(cl$par[2])
    K <- sdgp_cl^2 * exp(-outer(sort(unique(x)), sort(unique(x)), "-")^2 /
                           (2 * exp(2 * cl$par[3])))
    diag(K) <- diag(K) + 1e-12
    ch <- !inherits(try(chol(K), silent = TRUE), "try-error")
    bk <- Filter(function(k) identical(k$covstruct, "gp"), fr$frame$re_blocks)[[1]]
    cat(sprintf(paste0("NCP seed %d n %d A %g s %.0e | closed brms ll %.4f sdgp %.4g | ",
                       "centered 1e-6: ll %.4f sdgp %.4g conv %d %.1fs | ncp 1e-12: %s | ",
                       "chol(K + 1e-12) at brms max %s\n"),
                seed, n, A, s, -cl$objective, sdgp_cl, -fr$opt$objective,
                exp(fr$estimates$theta[bk$theta_idx][[1]]), fr$opt$convergence, t_c,
                r, ch))
  }
}
