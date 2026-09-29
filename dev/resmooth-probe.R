# Lane wt-resmooth. Which smooth blocks re_formula = NA keeps, in sample
# and on newdata, and what simulate(NA) redraws. Same data and seed as
# dev/simnewdata-review/rv-smooth.R so the numbers are comparable.
#   RESMOOTH_LIB=base Rscript dev/resmooth-probe.R > dev/resmooth-before.txt
#   Rscript dev/resmooth-probe.R > dev/resmooth-after.txt
lane_lib <- c("C:/Users/adf44/source/r/wt-resmooth-lib2",
              "C:/Users/adf44/source/r/wt-resmooth-lib")
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
.libPaths(c(if (!base) lane_lib,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)
d$ys <- sin(2 * pi * d$x) + rnorm(n, 0, exp(-1 + 1.5 * d$x))

nsim <- 1500
msg <- function(e) substr(conditionMessage(e), 1, 120)
tryv <- function(expr) tryCatch(expr, error = function(e) e)

probe <- function(label, f) {
  fit <- tryCatch(suppressWarnings(frm(f, data = d)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-22s FIT ERROR: %s\n", label, msg(fit))); return(invisible())
  }
  bl <- fit$frame$re_blocks
  cs <- vapply(bl, `[[`, "", "covstruct")
  b1 <- vapply(bl, function(b) b$b_idx[1], 0)
  kept <- function(nd) {
    k <- integer(0)
    for (lp in fit$frame$linpreds) {
      ed <- tryv(frmtmb:::lp_eta_design(fit, lp, nd, FALSE, FALSE))
      if (inherits(ed, "error")) return(paste0("ERR:", msg(ed)))
      k <- c(k, match(vapply(ed$sm_blocks, function(b) b$b_idx[1], 0), b1))
      # the newdata route reports its kept smooths through sm_parts
      if (!is.null(nd) && length(ed$sm_parts)) {
        k <- c(k, match(vapply(ed$sm_parts, function(p) p$bk$b_idx[1], 0), b1))
      }
    }
    paste(sort(unique(k)), collapse = ",")
  }
  ins <- kept(NULL)
  nds <- kept(d)
  redraw <- paste(frmtmb:::sim_re_plan(fit, NA)$blocks, collapse = ",")
  fna <- tryv(fitted(fit, re_formula = NA)[, "Estimate"])
  fnl <- tryv(fitted(fit, re_formula = NULL)[, "Estimate"])
  dnn <- if (inherits(fna, "error") || inherits(fnl, "error")) NA_real_ else
    max(abs(fna - fnl))
  # frm_linpred(), not predict(): predict.frmtmb_fit() SIMULATES the
  # predictive distribution, so its Estimate carries Monte Carlo noise
  # (0.089 on y ~ x, dev/resmooth-ndcheck.R) and cannot test a design
  pnd <- tryv(as.vector(frm_linpred(fit, newdata = d, re_formula = NA)))
  pin <- tryv(as.vector(frm_linpred(fit, re_formula = NA)))
  dnd <- if (inherits(pnd, "error") || inherits(pin, "error")) NA_real_ else
    max(abs(pnd - pin))
  s <- tryv(as.matrix(simulate(fit, nsim = nsim, seed = 2, re_formula = NA)))
  sgv <- as.vector(frm_linpred(fit, dpar = "sigma", type = "response",
                               re_formula = NA))
  if (length(sgv) == 1L) sgv <- rep(sgv, n)
  rsd <- if (inherits(s, "error")) NA_real_ else
    stats::median(apply(s, 1, sd) / sgv)
  # pp_check's own default: sd of the pooled draws over sd of the data
  pp <- tryv(as.matrix(simulate(fit, nsim = 200, seed = 3, re_formula = NA)))
  ppr <- if (inherits(pp, "error")) NA_real_ else
    sd(as.vector(pp)) / sd(fit$frame$y[[1L]])
  cat(sprintf(paste0("%-22s blocks %-28s | in-sample keeps %-7s | ",
                     "newdata keeps %-7s | sim redraws %-5s | ",
                     "max|NA-NULL| %9.3g | max|nd-ins| %9.3g | ",
                     "rowsd/sigma %.3f | sd(yrep)/sd(y) %.3f\n"),
              label, paste(seq_along(cs), cs, sep = ":", collapse = ","),
              ins, nds, redraw, dnn, dnd, rsd, ppr))
  if (inherits(pnd, "error")) cat("    predict(newdata) ERROR:",
                                  msg(pnd), "\n")
  if (inherits(s, "error")) cat("    simulate ERROR:", msg(s), "\n")
  invisible(fit)
}

probe("s(x)", bf(y ~ s(x)))
probe("s(g, bs = re)", bf(y ~ s(x) + s(g, bs = "re")))
probe("s(x, g, bs = fs)", bf(y ~ s(x, g, bs = "fs", k = 5)))
probe("t2(x, g, cr+re)", bf(y ~ t2(x, g, bs = c("cr", "re"))))
probe("s(x) + (1 | g)", bf(y ~ s(x) + (1 | g)))
probe("sigma ~ s(x)", bf(ys ~ s(x), sigma ~ s(x)))
probe("t2(x, z)", bf(y ~ t2(x, z)))
probe("gp(x)", bf(y ~ gp(x)))

cat("\n-- newdata with a NEW level of g, re_formula = NA --\n")
nd2 <- d[1:5, ]
nd2$g <- factor("99", levels = c(levels(d$g), "99"))
for (lab in c("s(g, bs = re)", "s(x, g, bs = fs)")) {
  f <- if (lab == "s(g, bs = re)") bf(y ~ s(x) + s(g, bs = "re")) else
    bf(y ~ s(x, g, bs = "fs", k = 5))
  fit <- suppressWarnings(frm(f, data = d))
  for (anl in c(FALSE, TRUE)) {
    p <- tryv(predict(fit, newdata = nd2, re_formula = NA,
                      allow_new_levels = anl))
    cat(sprintf("%-18s allow_new_levels=%-5s %s\n", lab, anl,
                if (inherits(p, "error")) paste("ERROR:", msg(p)) else
                  paste("OK:", paste(sprintf("%.4f", p[, "Estimate"]),
                                     collapse = " "))))
  }
  # and with the grouping column absent altogether
  nd3 <- nd2[, setdiff(names(nd2), "g"), drop = FALSE]
  p <- tryv(predict(fit, newdata = nd3, re_formula = NA))
  cat(sprintf("%-18s no g column        %s\n", lab,
              if (inherits(p, "error")) paste("ERROR:", msg(p)) else "OK"))
}

cat("\n-- conditional_effects(re_formula = NA) --\n")
for (lab in c("s(x, g, bs = fs)", "s(g, bs = re)")) {
  f <- if (lab == "s(g, bs = re)") bf(y ~ s(x) + s(g, bs = "re")) else
    bf(y ~ s(x, g, bs = "fs", k = 5))
  fit <- suppressWarnings(frm(f, data = d))
  ce <- tryv(conditional_effects(fit, effects = "x", resolution = 8))
  cat(sprintf("%-18s %s\n", lab,
              if (inherits(ce, "error")) paste("ERROR:", msg(ce)) else
                paste("OK range",
                      paste(sprintf("%.4f", range(ce[[1]]$estimate__)),
                            collapse = " "))))
}

cat("\n-- partial re_formula beside a factor smooth --\n")
fit <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | f)),
                            data = d))
p <- tryv(predict(fit, re_formula = ~ (1 | f)))
cat("partial:", if (inherits(p, "error")) paste("ERROR:", msg(p)) else "OK",
    "\n")
p <- tryv(predict(fit, re_formula = ~ (1 | nosuch)))
cat("nosuch :", if (inherits(p, "error")) paste("ERROR:", msg(p)) else "OK",
    "\n")
