# Where does frmtmb's nugget reach the fit? It is white noise of sd
# 1e-3 * sdgp on the latent field, so a gaussian fit can only tell it
# from sigma once sigma nears that size. Closed-form marginal maxima at
# each nugget, against frm(), over shrinking noise. Seed 20261102.
#
#   Rscript dev/gpgap-lownoise.R
source("dev/gpgap-lib.R")
E$local_run <- function() {
  set.seed(20261102)
  n <- 60
  x <- round(runif(n, 0, 6), 3)
  u <- sort(unique(x))
  D2 <- outer(u, u, "-")^2
  J <- match(x, u)
  e0 <- rnorm(n)
  for (s in c(1e-1, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)) {
    y <- sin(2 * x) + s * e0
    nll <- function(p, kind) {
      s2 <- exp(2 * p[2])
      K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
      diag(K) <- diag(K) + if (kind == "rel") 1e-6 * s2 else 1e-12
      V <- K[J, J]
      diag(V) <- diag(V) + exp(2 * p[4])
      R <- tryCatch(chol(V), error = function(e) NULL)
      if (is.null(R)) return(1e10)
      r <- backsolve(R, y - p[1], transpose = TRUE)
      0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
    }
    fit1 <- function(kind) {
      best <- NULL
      for (st in list(c(0, 0, 0, log(s)), c(0, 0, -0.5, log(s)),
                      c(0, 0.5, 0.5, log(s) - 1))) {
        o <- stats::nlminb(st, nll, kind = kind,
                           lower = c(-5, -5, -5, -25), upper = c(5, 5, 5, 3),
                           control = list(iter.max = 500, eval.max = 800))
        if (is.null(best) || o$objective < best$objective) best <- o
      }
      best
    }
    a <- fit1("rel")
    b <- fit1("abs")
    d <- data.frame(x = x, y = y)
    fr <- suppressWarnings(frmtmb::frm(frmtmb::bf(y ~ gp(x)), data = d))
    gp <- Filter(function(k) identical(k$covstruct, "gp"),
                 fr$frame$re_blocks)[[1]]
    sg <- as.numeric(stats::sigma(fr))
    cat(sprintf(paste0("GPLOW true sigma %.1e | closed form frm nugget: sigma ",
                       "%.4e sdgp %.4f rho %.4f ll %.4f | brms nugget: sigma ",
                       "%.4e sdgp %.4f rho %.4f ll %.4f | brms model at frm ",
                       "max loses %.3e | frm(): sigma %.4e sdgp %.4f conv %d\n"),
                s, exp(a$par[4]), exp(a$par[2]), exp(a$par[3]), -a$objective,
                exp(b$par[4]), exp(b$par[2]), exp(b$par[3]), -b$objective,
                nll(a$par, "abs") - b$objective, sg,
                exp(fr$estimates$theta[gp$theta_idx][[1]]),
                fr$opt$convergence))
    # the non-centered objective at brms's nugget, frm()'s own starts
    ns <- asNamespace("frmtmb")
    orig <- get("gp_corr", envir = ns)
    unlockBinding("gp_corr", ns)
    assign("gp_corr", function(theta, blk) {
      Q <- blk[["aux_D2"]][[1]] / (2 * exp(2 * theta[2]))
      exp(-Q) + diag(1, nrow(Q)) * (1e-12 * exp(-2 * theta[1]))
    }, envir = ns)
    frn <- fr$frame
    frn$ncp_blocks <- 1L
    fn <- frmtmb::build_objective(frn)
    tpl <- frn$par_template
    tpl$b[] <- 0
    r <- tryCatch({
      o <- RTMB::MakeADFun(fn, tpl, random = "b", map = frn$map, silent = TRUE)
      op <- suppressWarnings(stats::nlminb(o$par, o$fn, o$gr,
        control = list(iter.max = 500, eval.max = 800)))
      th <- op$par[names(op$par) == "theta"]
      sprintf("conv %d (%s) gr %.1e ll %.4f sdgp %.4f rho %.4f", op$convergence,
              substr(op$message, 1, 20), max(abs(o$gr(op$par))), -op$objective,
              exp(th[[1]]), exp(th[[2]]))
    }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 80)))
    assign("gp_corr", orig, envir = ns)
    lockBinding("gp_corr", ns)
    cat("  ncp at brms nugget:", r, "\n")
    # the Hilbert-space form has no nugget
    frk <- suppressWarnings(frmtmb::frm(frmtmb::bf(y ~ gp(x, k = 30)),
                                        data = d))
    hk <- Filter(function(k) identical(k$covstruct, "hsgp"),
                 frk$frame$re_blocks)[[1]]
    cat(sprintf("  gp(x, k = 30): sigma %.4e sdgp %.4f ll %.4f conv %d\n",
                as.numeric(stats::sigma(frk)),
                exp(frk$estimates$theta[hk$theta_idx][[1]]),
                -frk$opt$objective, frk$opt$convergence))
  }
}
environment(E$local_run) <- E
E$local_run()
