# Lane optstart: which change moves two standard-error tests on the lane
# build: test-se-tier.R's curved ridge (a nonlinear fit, whose start the
# lane does not change) and test-se-check.R's gr(g, by = f) fixture.
# Arms: the lane build as it is; with fd_hessian()'s step check off
# (tau_band = Inf, 0.70.0's Hessian); and for the gr() fixture also from
# 0.70.0's start, given as `start`. Per arm: the lost SEs and the
# warnings.
#   Rscript dev/optstart-sediag.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
fd0 <- get("fd_hessian", ns)
band <- function(on) {
  f <- if (on) fd0 else function(par, gr, tau_band = Inf, tau_step = 1e-4) {
    fd0(par, gr, tau_band = Inf, tau_step = tau_step)
  }
  environment(f) <- ns
  utils::assignInNamespace("fd_hessian", f, "frmtmb")
}
run <- function(label, expr) {
  w <- character()
  f <- withCallingHandlers(suppressMessages(expr), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf("%-28s code %d | lost %s | warnings: %s\n", label,
              f$opt$convergence, paste(names(lost), collapse = ","),
              substr(paste(w, collapse = " | "), 1, 160)))
  invisible(f)
}
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
fo <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
for (on in c(TRUE, FALSE)) {
  band(on)
  run(paste("ridge, band", on), frm(fo, data = dn))
}
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
uu <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + uu[dd$g, 1] +
                       uu[dd$g, 2] * dd$x, 1)
fg <- bf(y ~ x + (1 + x | gr(g, by = f)))
old <- list(beta = c(mean(dd$y), 0), betad = log(stats::sd(dd$y)))
for (on in c(TRUE, FALSE)) {
  band(on)
  run(paste("gr by, band", on), frm(fg, family = gaussian(), data = dd))
  run(paste("gr by old start, band", on),
      frm(fg, family = gaussian(), data = dd, start = old))
}
