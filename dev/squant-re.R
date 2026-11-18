# Lane squant, item 4(d): quantile_elf() with (1 | g) in mu. On the
# designs of backlog item 5449 (dev/reviews/wcheck2-review-al.R, seeds
# 20261114 to 20261123), for each smoothing: the Laplace fit (code,
# warnings, logLik, log sd), the Laplace objective profiled along
# log sd (every other parameter re-optimized), and the quadrature fit.
#
# Run: Rscript dev/squant-re.R <design> <quantile>
#   design: g4, g10, g2, al4 (as in wcheck2-review-al.R)
# Output on stdout; the driver writes dev/squant-out/re-<design>-<q>.txt
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
a <- commandArgs(TRUE)
design <- a[1]
tau <- as.numeric(a[2])
hs <- c(0.1, 0.25, 0.5, 1)
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")

mk_data <- function(seed, design) {
  set.seed(seed)
  gs <- switch(design, g4 = 4, al4 = 4, g10 = 10, g2 = 2)
  n <- if (design == "g10") 100 else 80
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / gs), each = gs)))
  u <- rnorm(n / gs, 0, 0.7)[as.integer(d$g)]
  e <- if (design == "al4") {
    p <- 0.3; v <- rexp(n); w <- rexp(n)
    0.5 * (v / p - w / (1 - p))
  } else rnorm(n)
  d$y <- 1 + 0.5 * d$x + u + e
  d
}

quiet <- function(expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(gsub("\\s+", " ", conditionMessage(w)), 1, 70))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  list(value = r, warnings = ws)
}

# the Laplace objective along log sd, every other outer parameter
# re-optimized from the fit; returns the log-likelihood at each offset
profile_sd <- function(fit, offsets) {
  obj <- fit$obj
  par <- fit$opt$par
  it <- which(names(par) == "theta")[1L]
  vapply(offsets, function(k) {
    th <- par[it] + k
    f <- function(p) {
      full <- par
      full[-it] <- p
      full[it] <- th
      v <- obj$fn(full)
      if (!is.finite(v)) 1e10 else v
    }
    o <- nlminb(par[-it], f, control = list(iter.max = 300, eval.max = 600))
    -o$objective
  }, 0)
}

offsets <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8)
fo <- bf(y ~ x + (1 | g), quantile = tau)
rows <- list()
for (s in 1:10) {
  d <- mk_data(20261113 + s, design)
  for (h in hs) {
    fam <- quantile_elf(smoothing = h)
    lf <- quiet(frm(fo, family = fam, data = d))
    if (inherits(lf$value, "error")) {
      cat(sprintf("seed %d h %.2f: Laplace ERROR %s\n", s, h,
                  conditionMessage(lf$value)))
      next
    }
    fit <- lf$value
    ll <- as.numeric(logLik(fit))
    lsd <- fit$opt$par[names(fit$opt$par) == "theta"][1L]
    pr <- profile_sd(fit, offsets)
    # an interior maximum: the fit's value is the profile's largest, to
    # the optimizer's tolerance, and the ends fall below it
    top <- max(pr)
    interior <- pr[offsets == 0] >= top - 1e-4 &&
      pr[1L] < top - 1e-3 && pr[length(pr)] < top - 1e-3
    qf <- quiet(frm(fo, family = fam, data = d, quadrature = TRUE))
    if (inherits(qf$value, "error")) {
      llq <- NA
      lsdq <- NA
      qnote <- substr(conditionMessage(qf$value), 1, 70)
    } else {
      llq <- as.numeric(logLik(qf$value))
      lsdq <- qf$value$opt$par[names(qf$value$opt$par) == "theta"][1L]
      qnote <- paste(qf$warnings, collapse = " | ")
    }
    cat(sprintf(paste0("seed %d h %.2f: conv %d warn %d logLik %.4f log sd",
                       " %.4f | profile max at offset %+.1f, ends %.3f %.3f",
                       " below, interior %s | quad logLik %.4f log sd %.4f",
                       " conv %s %s%s\n"),
                s, h, fit$opt$convergence, length(lf$warnings), ll, lsd,
                offsets[which.max(pr)], top - pr[1L], top - pr[length(pr)],
                interior, llq, lsdq,
                if (is.list(qf$value) && !inherits(qf$value, "error"))
                  qf$value$opt$convergence else "NA",
                qnote,
                if (length(lf$warnings)) paste0(" [",
                  paste(unique(lf$warnings), collapse = " | "), "]") else ""))
    cat("  profile:", paste(sprintf("%.4f", pr - top), collapse = " "), "\n")
    rows[[length(rows) + 1L]] <- data.frame(
      seed = s, h = h, conv = fit$opt$convergence,
      warn = length(lf$warnings), interior = interior, ll = ll, lsd = lsd,
      llq = llq, lsdq = lsdq)
  }
}
r <- do.call(rbind, rows)
for (h in hs) {
  x <- r[r$h == h, ]
  cat(sprintf(paste0("SUMMARY %s q %.2f h %.2f: %d fits, conv 0 %d, no warning",
                     " %d, interior profile maximum %d; quadrature fits %d;",
                     " |log sd - quad| median %.4f max %.4f; logLik - quad",
                     " median %.4f range %.4f %.4f\n"),
              design, tau, h, nrow(x), sum(x$conv == 0), sum(x$warn == 0),
              sum(x$interior), sum(is.finite(x$llq)),
              median(abs(x$lsd - x$lsdq), na.rm = TRUE),
              max(abs(x$lsd - x$lsdq), na.rm = TRUE),
              median(x$ll - x$llq, na.rm = TRUE),
              min(x$ll - x$llq, na.rm = TRUE),
              max(x$ll - x$llq, na.rm = TRUE)))
}
cat("DONE\n")
