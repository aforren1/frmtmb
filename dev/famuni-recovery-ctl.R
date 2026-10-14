# The control for dev/famuni-recovery.R: the same design, seeds and
# random intercept under gaussian(), so that a standard-error shortfall
# seen there can be read as the families' own or as the design's.
# Usage: Rscript dev/famuni-recovery-ctl.R gauss [R]
# Output: dev/famuni-out/recovery-ctl-gauss-<R>.txt
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
which_fam <- commandArgs(TRUE)[1]
R <- if (length(commandArgs(TRUE)) > 1L) as.integer(commandArgs(TRUE)[2]) else 40L
n <- 400L
G <- 20L
sim1 <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(G, n / G))
  b <- rnorm(G, 0, 0.3)[d$g]
  if (which_fam == "frechet") {
    nu <- 1 + exp(log(3) + 0.2 * d$z)
    d$y <- exp(1 + 0.3 * d$x + b) / gamma(1 - 1 / nu) * rexp(n)^(-1 / nu)
    list(d = d, f = bf(y ~ x + (1 | g), nu ~ z), fam = frechet(),
         truth = c(Intercept = 1, x = 0.3, nu_Intercept = log(3),
                   nu_z = 0.2, sd = 0.3))
  } else if (which_fam == "gev") {
    xi <- 0.15
    sg <- exp(-0.3 + 0.2 * d$z)
    d$y <- 1 + 0.5 * d$x + b + sg * (rexp(n)^(-xi) - 1) / xi
    list(d = d, f = bf(y ~ x + (1 | g), sigma ~ z), fam = gen_extreme_value(),
         truth = c(Intercept = 1, x = 0.5, sigma_Intercept = -0.3,
                   sigma_z = 0.2, xi = xi, sd = 0.3))
  } else {
    # the control: the same design and random intercept under gaussian(),
    # whose standard errors are long established
    d$y <- 1 + 0.5 * d$x + b + exp(-0.3 + 0.2 * d$z) * rnorm(n)
    list(d = d, f = bf(y ~ x + (1 | g), sigma ~ z), fam = gaussian(),
         truth = c(Intercept = 1, x = 0.5, sigma_Intercept = -0.3,
                   sigma_z = 0.2, sd = 0.3))
  }
}
rows <- list()
for (s in seq_len(R)) {
  sm <- sim1(s)
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frm(sm$f, family = sm$fam, data = sm$d),
             error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (inherits(fit, "error")) {
    cat("seed", s, "ERROR", conditionMessage(fit), "\n")
    next
  }
  cf <- fixef(fit)
  est <- setNames(cf[, "Estimate"], rownames(cf))
  se <- setNames(cf[, "Est.Error"], rownames(cf))
  if ("xi" %in% names(sm$truth)) {
    lp <- frm_linpred(fit, dpar = "xi", type = "response", se.fit = TRUE)
    est["xi"] <- lp$fit[1]
    se["xi"] <- lp$se.fit[1]
  }
  vc <- VarCorr(fit)
  est["sd"] <- vc[["g"]][["sd"]][1, "Estimate"]
  se["sd"] <- vc[["g"]][["sd"]][1, "Est.Error"]
  conv <- fit$opt$convergence %||% NA
  pdh <- all(is.finite(se))
  nse <- sum(!is.finite(se[names(sm$truth) != "sd"]))
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, par = names(sm$truth), truth = unname(sm$truth),
    est = unname(est[names(sm$truth)]), se = unname(se[names(sm$truth)]),
    conv = conv, pdHess = pdh, nwarn = length(warns),
    warn = paste(unique(warns), collapse = " | "))
}
out <- do.call(rbind, rows)
cat("family", which_fam, "seeds with a fit:", length(unique(out$seed)), "of",
    R, "\n")
cat("convergence code 0:", sum(tapply(out$conv, out$seed, `[`, 1) == 0),
    "; all SEs finite:", sum(tapply(out$pdHess, out$seed, `[`, 1)),
    "; fits with any warning:", sum(tapply(out$nwarn, out$seed, `[`, 1) > 0),
    "\n")
w <- unique(out$warn[nzchar(out$warn)])
if (length(w)) cat("warnings seen:\n", paste(" ", w, collapse = "\n"), "\n")
cat(sprintf("%-16s %7s %9s %9s %9s %9s %9s\n", "parameter", "truth",
            "mean est", "bias", "emp SD", "mean SE", "cover95"))
for (p in unique(out$par)) {
  o <- out[out$par == p, ]
  cov <- if (all(is.na(o$se))) NA else
    mean(abs(o$est - o$truth) <= 1.96 * o$se, na.rm = TRUE)
  cat(sprintf("%-16s %7.3f %9.4f %9.4f %9.4f %9.4f %9s\n", p, o$truth[1],
              mean(o$est), mean(o$est) - o$truth[1], sd(o$est),
              mean(o$se, na.rm = TRUE),
              if (is.na(cov)) "-" else sprintf("%.3f", cov)))
}
