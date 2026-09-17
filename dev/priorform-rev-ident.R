# Reviewer, lane wt-priorform: are the frmtmb-only designs the lane now
# refuses identified? Each is simulated from known values and fitted on
# the BASE build, which accepts them. Seeds 1..4 per design.
#   Rscript dev/priorform-rev-ident.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
ng <- 150; nt <- 4
sim <- function(seed, kind) {
  set.seed(seed)
  d <- expand.grid(g = factor(1:ng), time = factor(1:nt))
  d$pos <- num_factor(as.integer(d$time))
  d$f <- d$time
  S <- switch(kind,
    ar1_diag = 1.0 * 0.7^abs(outer(1:nt, 1:nt, "-")) + diag(c(.3, .6, .9, 1.2)),
    rr_diag = tcrossprod(c(1, .8, .6, 1.2)) + diag(c(.3, .5, .4, .6)),
    equalto_us = { V <- matrix(.2, nt, nt) + diag(.3, nt)
      U <- matrix(.3, nt, nt) + diag(.5, nt); V + U },
    exp_diag = exp(-abs(outer(1:nt, 1:nt, "-")) / 2) + diag(.4, nt),
    us_diag = matrix(.3, nt, nt) + diag(.8, nt),
    gr_plain = NULL)
  if (kind == "gr_plain") {
    lev <- levels(d$g)
    A <- diag(ng)
    for (i in seq(1, ng - 1, 2)) A[i, i + 1] <- A[i + 1, i] <- 0.5
    dimnames(A) <- list(lev, lev)
    u <- drop(t(chol(A)) %*% rnorm(ng)) * 0.8
    pe <- rnorm(ng) * 0.6
    d$y <- u[d$g] + pe[d$g] + rnorm(nrow(d), 0, 0.5)
    return(list(d = d, A = A))
  }
  b <- t(chol(S) %*% matrix(rnorm(ng * nt), nt, ng))
  d$y <- b[cbind(as.integer(d$g), as.integer(d$time))] + rnorm(nrow(d), 0, 0.3)
  list(d = d, V = matrix(.2, nt, nt) + diag(.3, nt))
}
forms <- list(
  ar1_diag = "bf(y ~ 1 + ar1(0 + time | g) + diag(0 + time | g))",
  rr_diag = "bf(y ~ 1 + rr(0 + time | g, d = 1) + diag(0 + time | g))",
  equalto_us = "bf(y ~ 1 + equalto(0 + f | g, V) + (0 + f | g))",
  exp_diag = "bf(y ~ 1 + exp(0 + pos | g) + diag(0 + pos | g))",
  us_diag = "bf(y ~ 1 + us(0 + time | g) + diag(0 + time | g))",
  gr_plain = "bf(y ~ 1 + (1 | gr(g, cov = A)) + (1 | g))")
for (k in names(forms)) {
  for (seed in 1:4) {
    s <- sim(seed, k)
    warns <- character(0)
    r <- tryCatch(withCallingHandlers({
      fit <- suppressMessages(frm(eval(parse(text = forms[[k]])) + gaussian(),
                                  data = s$d,
                                  data2 = list(A = s$A, V = s$V)))
      th <- fit$opt$par[names(fit$opt$par) == "theta"]
      sdr <- frmtmb:::sdr_of(fit)
      cf <- sdr$cov.fixed
      sprintf("conv=%d pdHess=%s ll=%.4f theta=[%s] max_theta_se=%.3g",
              fit$opt$convergence, isTRUE(sdr$pdHess),
              as.numeric(logLik(fit)),
              paste(sprintf("%.2f", th), collapse = " "),
              max(sqrt(diag(cf))[rownames(cf) == "theta"]))
    }, warning = function(w) {
      warns <<- c(warns, substr(conditionMessage(w), 1, 50))
      invokeRestart("muffleWarning")
    }), error = function(e) paste("ERROR:", conditionMessage(e)))
    if (length(warns)) {
      r <- paste(r, "| warned:", paste(unique(warns), collapse = "; "))
    }
    cat(k, seed, r, "\n")
  }
}
