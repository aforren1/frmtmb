# Review of lane squant, item 5: with (1 | g), the filed MEDIUM (a fit at
# code 0 below a higher point of its own Laplace profile in log sd) on
# 5449's design g4, quantile 0.9, seed 10, smoothing 0.25 and 0.1; a
# refit from the profile's best point; the far ends of the profile;
# and the sandwich and quadrature routes the docs recommend.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-re.R > dev/reviews/squant-review-re.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
mk_data <- function(seed, gs = 4, n = 80) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / gs), each = gs)))
  u <- rnorm(n / gs, 0, 0.7)[as.integer(d$g)]
  d$y <- 1 + 0.5 * d$x + u + rnorm(n)
  d
}
quiet <- function(expr) {
  ws <- character()
  v <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(gsub("\\s+", " ", conditionMessage(w)), 1, 120))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ws <<- c(ws, paste("MESSAGE", substr(conditionMessage(m), 1, 120)))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  list(v = v, w = ws)
}
profile_sd <- function(fit, offsets) {
  obj <- fit$obj
  par <- fit$opt$par
  it <- which(names(par) == "theta")[1L]
  lapply(offsets, function(k) {
    f <- function(p) {
      full <- par
      full[-it] <- p
      full[it] <- par[it] + k
      v <- obj$fn(full)
      if (!is.finite(v)) 1e10 else v
    }
    o <- nlminb(par[-it], f, control = list(iter.max = 500, eval.max = 1000))
    full <- par
    full[-it] <- o$par
    full[it] <- par[it] + k
    list(ll = -o$objective, par = full)
  })
}
offsets <- c(-6, -3, -2, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 4, 8, 12)
fo <- bf(y ~ x + (1 | g), quantile = 0.9)
d <- mk_data(20261113 + 10)
for (h in c(0.1, 0.25, 1)) {
  r <- quiet(frm(fo, family = quantile_elf(smoothing = h), data = d))
  fit <- r$v
  pr <- profile_sd(fit, offsets)
  ll <- vapply(pr, `[[`, 0, "ll")
  top <- which.max(ll)
  cat(sprintf(paste0("h %.2f: conv %d, %d conditions, logLik %.4f, log sd",
                     " %.4f | profile - fit at offsets: %s\n"),
              h, fit$opt$convergence, length(r$w), -fit$opt$objective,
              fit$opt$par[names(fit$opt$par) == "theta"][1],
              paste(sprintf("%+g:%.3f", offsets, ll + fit$opt$objective),
                    collapse = " ")))
  if (length(r$w)) cat("   conditions:", paste(r$w, collapse = " | "), "\n")
  # a fit started at the profile's best point: does frm() go there?
  st <- pr[[top]]$par
  o2 <- nlminb(st, fit$obj$fn, fit$obj$gr,
               control = list(iter.max = 1000, eval.max = 2000))
  cat(sprintf(paste0("   nlminb from the profile's best point (offset %+g):",
                     " logLik %.4f, log sd %.4f, code %d (fit's %.4f)\n"),
              offsets[top], -o2$objective,
              o2$par[names(o2$par) == "theta"][1], o2$convergence,
              -fit$opt$objective))
  q <- quiet(frm(fo, family = quantile_elf(smoothing = h), data = d,
                 quadrature = TRUE))
  if (inherits(q$v, "error")) {
    cat("   quadrature: ERROR", conditionMessage(q$v), "\n")
  } else {
    cat(sprintf("   quadrature: conv %d, logLik %.4f, log sd %.4f, %d conditions\n",
                q$v$opt$convergence, as.numeric(logLik(q$v)),
                q$v$opt$par[names(q$v$opt$par) == "theta"][1], length(q$w)))
    s <- quiet(vcov_cluster(q$v, ~ g, "CR0"))
    cat("   vcov_cluster on the quadrature fit:",
        if (inherits(s$v, "error")) paste("ERROR", substr(conditionMessage(s$v),
                                                          1, 200))
        else paste(sprintf("%.4f", sqrt(diag(s$v))), collapse = " "), "\n")
  }
  s <- quiet(vcov_cluster(fit, ~ g, "CR0"))
  cat("   vcov_cluster(~ g) on the Laplace fit:",
      if (inherits(s$v, "error")) paste("ERROR", substr(conditionMessage(s$v),
                                                        1, 200))
      else paste(sprintf("%.4f", sqrt(diag(s$v))), collapse = " "), "\n")
}

# (1 | obs) with unit weights: the exact marginal is bounded (ELF
# convolved with a normal); does the Laplace fit end at a maximum or run?
cat("\n(1 | obs), unit weights, smoothing 0.25 and 1, 3 seeds\n")
for (s in 1:3) {
  for (h in c(0.25, 1)) {
    set.seed(s)
    n <- 60
    dd <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
    dd$y <- 1 + 0.5 * dd$x + rnorm(n)
    r <- quiet(frm(bf(y ~ x + (1 | obs), quantile = 0.5),
                   family = quantile_elf(smoothing = h), data = dd))
    if (inherits(r$v, "error")) {
      cat(sprintf("  seed %d h %.2f: ERROR %s\n", s, h,
                  substr(conditionMessage(r$v), 1, 160)))
      next
    }
    f <- r$v
    par <- f$opt$par
    is <- grep("sigma|betad", names(par))[1]
    prof <- vapply(c(-2, -5, -10), function(k) {
      o <- nlminb(par[-is], function(p) {
        full <- par
        full[-is] <- p
        full[is] <- par[is] + k
        v <- f$obj$fn(full)
        if (!is.finite(v)) 1e10 else v
      })
      -o$objective
    }, 0)
    cat(sprintf(paste0("  seed %d h %.2f: conv %d, %d conditions, logLik %.3f,",
                       " sigma par %s %.3f | profile at -2 -5 -10: %s %s\n"),
                s, h, f$opt$convergence, length(r$w), -f$opt$objective,
                names(par)[is], par[is],
                paste(sprintf("%.3f", prof), collapse = " "),
                paste(r$w, collapse = " | ")))
  }
}
cat("DONE\n")
