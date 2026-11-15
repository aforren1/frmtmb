# Review of lane wcheck2: how large and how common is the inner-search
# error of asym_laplace with a grouped effect, no weights, no censoring?
# For each design and seed: frm()'s fit (code, warnings, logLik), the
# Laplace approximation with each group's exact inner mode at the fit's
# own parameters, the optimum of that exact-mode Laplace approximation
# (outer parameters re-optimized here), and the quadrature fit.
#   Rscript dev/reviews/wcheck2-review-al.R <base|lane> <design>
# design: g4 (80 rows, groups of 4), g10 (100 rows, groups of 10),
#         g2 (80 rows, groups of 2), al4 (g4 with ALD noise),
#         huber4 (huber, groups of 4: inner gradient at TMB's modes)
a <- commandArgs(TRUE)
lib_arg <- a[1]
design <- a[2]
source("dev/wcheck2-common.R")

mk_data <- function(seed, design) {
  set.seed(seed)
  gs <- switch(design, g4 = 4, al4 = 4, huber4 = 4, g10 = 10, g2 = 2)
  n <- if (design == "g10") 100 else 80
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / gs), each = gs)))
  u <- rnorm(n / gs, 0, 0.7)[as.integer(d$g)]
  e <- if (design == "al4") {
    # asymmetric Laplace noise at quantile 0.3, scale 0.5
    p <- 0.3; v <- rexp(n); w <- rexp(n)
    0.5 * (v / p - w / (1 - p))
  } else rnorm(n)
  d$y <- 1 + 0.5 * d$x + u + e
  d
}

ald_lp <- function(y, mu, s, p) {
  u <- (y - mu) / s
  log(p * (1 - p) / s) - u * (p - (u < 0))
}

# exact-mode Laplace of the model at outer parameters
lap_exact <- function(par, d, nm) {
  beta <- par[nm == "beta"]
  bd <- par[nm == "betad"]
  s <- exp(bd[1L]); p <- plogis(bd[2L])
  sd <- exp(par[nm == "theta"][1L])
  eta <- beta[1L] + beta[2L] * d$x
  grp <- as.integer(d$g)
  G <- max(grp)
  tot <- 0
  modes <- numeric(G)
  for (k in seq_len(G)) {
    i <- which(grp == k)
    f <- function(b) sum(ald_lp(d$y[i], eta[i] + b, s, p)) - b^2 / (2 * sd^2)
    # the maximum of a concave piecewise quadratic sits at a kink or in
    # a smooth piece: optimize() over a wide bracket, then the kinks
    o <- optimize(f, c(-20 * sd - 20, 20 * sd + 20), maximum = TRUE,
                  tol = 1e-12)
    kinks <- d$y[i] - eta[i]
    fk <- vapply(kinks, f, 0)
    if (max(fk) > o$objective) {
      o$objective <- max(fk); o$maximum <- kinks[which.max(fk)]
    }
    modes[k] <- o$maximum
    # data curvature is 0 away from the kinks: H = 1 / sd^2, and the
    # Laplace term is f(b*) - log N's constant cancel: f + log(sd) +
    # 0.5 log(2 pi) - 0.5 log(1 / sd^2) - log(sd) - 0.5 log(2 pi)
    tot <- tot + o$objective
  }
  list(val = tot, modes = modes)
}

seeds <- 1:10
cat("design", design, "\n")
rows <- list()
for (sd0 in seeds) {
  d <- mk_data(20261113 + sd0, design)
  fam <- if (design == "huber4") huber() else asym_laplace()
  fo <- y ~ x + (1 | g)
  ws <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(frm(fo, family = fam, data = d),
    warning = function(w) {
      ws <<- c(ws, paste0(class(w)[1], ": ", substr(gsub("\\s+", " ",
                                                         conditionMessage(w)),
                                                    1, 60)))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  tfit <- proc.time()[["elapsed"]] - t0
  par <- fit$opt$par
  nm <- names(par)
  ll <- as.numeric(logLik(fit))
  ob <- q(frm(fo, family = fam, data = d, dry_run = "objective"))
  obj <- if (!is.null(ob$obj)) ob$obj else ob
  tmb_fresh <- -obj$fn(par)
  if (design == "huber4") {
    # inner gradient of the joint at TMB's own modes, a smooth density
    gfull <- obj$env$f(obj$env$last.par, order = 1)
    gin <- max(abs(gfull[obj$env$random]))
    cat(sprintf(paste0("seed %d: conv %d, warn %d, logLik %.5f, fresh TMB %.5f,",
                       " max |inner gradient| %.3g%s\n"),
                sd0, fit$opt$convergence, length(ws), ll, tmb_fresh, gin,
                if (length(ws)) paste0(" [", paste(ws, collapse = " | "), "]")
                else ""))
    next
  }
  le <- lap_exact(par, d, nm)
  # TMB's modes against the exact ones, at the fit's parameters
  bh <- obj$env$last.par[obj$env$random]
  gap_mode <- max(abs(bh - le$modes))
  # the optimum of the exact-mode Laplace approximation
  op <- nlminb(par, function(p) -lap_exact(p, d, nm)$val,
               control = list(iter.max = 300, eval.max = 600))
  op2 <- nlminb(op$par, function(p) -lap_exact(p, d, nm)$val,
                control = list(iter.max = 300, eval.max = 600))
  if (op2$objective < op$objective) op <- op2
  qerr <- ""
  fq <- tryCatch(q(frm(fo, family = fam, data = d, quadrature = TRUE)),
                 error = function(e) {
                   qerr <<- substr(gsub("\\s+", " ", conditionMessage(e)),
                                   1, 90)
                   NULL
                 })
  llq <- if (is.null(fq)) NA else as.numeric(logLik(fq))
  if (nzchar(qerr)) cat("  quadrature error:", qerr, "\n")
  # the exact-mode Laplace along log sd, the rest at the fit: it never
  # falls as sd grows (data curvature 0, so -log H cancels the prior's
  # normalizing constant)
  it <- which(nm == "theta")[1L]
  prof_sd <- vapply(c(0, 1, 2, 4, 8), function(k) {
    pp <- par; pp[it] <- par[it] + k; lap_exact(pp, d, nm)$val
  }, 0)
  cat("  exact-mode Laplace at log sd + 0, 1, 2, 4, 8:",
      paste(format(prof_sd, digits = 7), collapse = " "), "\n")
  # the fit's dispersion and sd against the exact-mode optimum's
  dpar <- par - op$par
  cat(sprintf(paste0("seed %d: conv %d, warn %d, logLik %.4f | fresh TMB at",
                     " par %.4f | exact-mode Laplace at par %.4f | its",
                     " optimum %.4f (gap %.4f) | quad logLik %.4f | max",
                     " |mode gap| %.3g | par diff beta %s, log sigma %.3f,",
                     " logit q %.3f, log sd %.3f | fit %.1fs%s\n"),
              sd0, fit$opt$convergence, length(ws), ll, tmb_fresh, le$val,
              -op$objective, -op$objective - ll, llq, gap_mode,
              paste(format(dpar[nm == "beta"], digits = 3), collapse = " "),
              dpar[nm == "betad"][1], dpar[nm == "betad"][2],
              dpar[nm == "theta"][1], tfit,
              if (length(ws)) paste0(" [", paste(unique(ws), collapse = " | "),
                                     "]") else ""))
  rows[[length(rows) + 1L]] <- data.frame(seed = sd0,
    conv = fit$opt$convergence, warn = length(ws), gap = -op$objective - ll,
    silent = fit$opt$convergence == 0 && !length(ws))
}
if (length(rows)) {
  r <- do.call(rbind, rows)
  cat(sprintf(paste0("SUMMARY %s: %d fits, conv 0 %d, no warning %d, silent",
                     " (conv 0, no warning) %d; logLik gap to the exact-mode",
                     " Laplace optimum: median %.3f, max %.3f; silent max",
                     " %.3f\n"),
              design, nrow(r), sum(r$conv == 0), sum(r$warn == 0),
              sum(r$silent), median(r$gap), max(r$gap),
              if (any(r$silent)) max(r$gap[r$silent]) else NA))
}
cat("DONE\n")
