# What does frmtmb's nugget do to the fit, and can frmtmb fit at brms's?
#
#   Rscript dev/gpgap-nugvar.R "<spec key>"
#
# Refits the pair (any prior dropped) with gp_corr() replaced in the
# namespace of the base build: relative nuggets 1e-6 (shipped), 1e-9,
# 1e-12 and brms's absolute 1e-12, on the centered objective frm() uses
# and, for the absolute one, on the non-centered objective (b = L z).
# Where the family lets the field be integrated in closed form (gaussian,
# lognormal; y1 of a multivariate fit) the same variants are maximized on
# y's marginal density, whose covariance is conditioned by sigma^2 and
# needs no inverse of the kernel: that is the truth for each nugget.
args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args)) args[1] else "lognormal x gp"
source("dev/gpgap-lib.R")
E$key <- key
local(envir = E, {
  ns <- asNamespace("frmtmb")
  corr_src <- function(kind, nu) {
    function(theta, blk) {
      Q <- 0
      for (j in seq_along(blk[["aux_D2"]])) {
        rho <- if (isTRUE(blk[["gp_iso"]])) exp(theta[2]) else exp(theta[1 + j])
        Q <- Q + blk[["aux_D2"]][[j]] / (2 * rho^2)
      }
      n <- nrow(blk[["aux_D2"]][[1]])
      I <- diag(1, n)
      if (identical(kind, "rel")) {
        exp(-Q) + I * nu
      } else {
        exp(-Q) + I * (nu * exp(-2 * theta[1]))
      }
    }
  }
  patch <- function(fn) {
    unlockBinding("gp_corr", ns)
    assign("gp_corr", fn, envir = ns)
    lockBinding("gp_corr", ns)
  }
  spec <- eq_specs[[key]]()
  spec$prior <- NULL
  orig <- get("gp_corr", envir = ns)

  gp_read <- function(fit) {
    bks <- Filter(function(b) identical(b$covstruct, "gp"),
                  fit$frame$re_blocks)
    th <- fit$estimates$theta[bks[[1]]$theta_idx]
    c(sdgp = exp(th[[1]]), rho = exp(th[[2]]))
  }
  one <- function(lab, fn) {
    patch(fn)
    t0 <- proc.time()[[3]]
    fit <- tryCatch(suppressWarnings(eq_frm_fit(spec)),
                    error = function(e) conditionMessage(e))
    dt <- proc.time()[[3]] - t0
    if (is.character(fit)) {
      cat(sprintf("  %-14s ERROR %s (%.1fs)\n", lab, substr(fit, 1, 70), dt))
      return(NULL)
    }
    g <- gp_read(fit)
    gr <- tryCatch(max(abs(fit$obj$gr(fit$opt$par))), error = function(e) NA)
    cat(sprintf(paste0("  %-14s conv %d (%s) gr %.2e, ll %.6f, sdgp %.5f, ",
                       "rho %.5f (%.1fs)\n"),
                lab, fit$opt$convergence, substr(fit$opt$message, 1, 24),
                gr, -fit$opt$objective, g[["sdgp"]], g[["rho"]], dt))
    fit
  }
  cat("GPVAR", key, "\n")
  f6 <- one("rel 1e-6", corr_src("rel", 1e-6))
  one("rel 1e-9", corr_src("rel", 1e-9))
  one("rel 1e-12", corr_src("rel", 1e-12))
  one("abs 1e-12", corr_src("abs", 1e-12))

  # the non-centered objective at brms's nugget, from frm()'s own starts
  patch(corr_src("abs", 1e-12))
  frn <- f6$frame
  frn$ncp_blocks <- which(vapply(frn$re_blocks,
                                 function(b) identical(b$covstruct, "gp"), NA))
  fn <- frmtmb::build_objective(frn)
  bk <- Filter(function(b) identical(b$covstruct, "gp"), frn$re_blocks)[[1]]
  ncp_run <- function(lab, tpl) {
    t0 <- proc.time()[[3]]
    r <- tryCatch({
      o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
      op <- suppressWarnings(stats::nlminb(o$par, o$fn, o$gr))
      th <- op$par[names(op$par) == "theta"][bk$theta_idx]
      sprintf("conv %d (%s) gr %.2e, ll %.6f, sdgp %.5f, rho %.5f",
              op$convergence, substr(op$message, 1, 24),
              max(abs(o$gr(op$par))), -op$objective, exp(th[[1]]), exp(th[[2]]))
    }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 70)))
    cat(sprintf("  %-22s %s (%.1fs)\n", lab, r, proc.time()[[3]] - t0))
  }
  tpl <- frn$par_template
  tpl$b[] <- 0
  ncp_run("abs 1e-12 ncp, default", tpl)
  tpl$theta <- f6$estimates$theta
  ncp_run("abs 1e-12 ncp, at fit", tpl)
  patch(orig)

  # Closed form
  d <- spec$data
  yn <- if (grepl("y1", spec$bf)) "y1" else "y"
  fam <- spec$family %||% "gaussian()"
  if (!grepl("student", fam) && !grepl("nl = TRUE", spec$bf)) {
    y <- d[[yn]]
    if (grepl("lognormal", fam)) y <- log(y)
    u <- sort(unique(d$x))
    J <- match(d$x, u)
    D2 <- outer(u, u, "-")^2
    nll <- function(p, kind, nu) {
      s2 <- exp(2 * p[2])
      K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
      diag(K) <- diag(K) + if (kind == "rel") nu * s2 else nu
      V <- K[J, J]
      diag(V) <- diag(V) + exp(2 * p[4])
      R <- chol(V)
      r <- backsolve(R, y - p[1], transpose = TRUE)
      0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
    }
    g6 <- gp_read(f6)
    s0 <- exp(f6$estimates$beta_sigma[[1]] %||% 0)
    p0 <- c(mean(y), log(g6[["sdgp"]]), log(g6[["rho"]]), log(stats::sd(y) / 2))
    best <- list()
    for (v in list(c("rel", 1e-6), c("abs", 1e-12))) {
      o <- stats::nlminb(p0, nll, kind = v[1], nu = as.numeric(v[2]),
                         control = list(iter.max = 500, rel.tol = 1e-12))
      best[[v[1]]] <- o
      cat(sprintf(paste0("  closed form %s %s: ll %.6f, sdgp %.5f, rho %.5f, ",
                         "sigma %.5f, conv %d\n"), v[1], v[2], -o$objective,
                  exp(o$par[2]), exp(o$par[3]), exp(o$par[4]), o$convergence))
    }
    cat(sprintf("  brms's model read at frmtmb's maximum: ll %.6f (loss %.3e)\n",
                -nll(best$rel$par, "abs", 1e-12),
                nll(best$rel$par, "abs", 1e-12) - best$abs$objective))
  }
})
