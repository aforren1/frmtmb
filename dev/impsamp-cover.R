# Is the Monte Carlo error the fit reports the error it makes? Over draw
# seeds 1:40, z = (logLik - exact at the fit's own estimates) / MCSE, with
# the exact value from the independent reference. A calibrated MCSE has
# sd(z) near 1. Two MCSEs are read from the same weights: the iid formula
# imp_ess() has used (a pair of antithetic draws counted as two draws)
# and one that treats each antithetic pair as the unit.
# Designs: p124 (seed 124), the hard design (seed 3), a poisson correlated
# slope (grid data, 40x5, seed 2). Log: dev/impsamp-out/cover.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
source("dev/impsamp-agh.R")
E <- new.env(parent = ns)
sys.source("dev/equiv-specs.R", envir = E)
d1 <- E$eq_dat("binomial", seed = 124)
set.seed(3)
ng <- 40L
g <- factor(rep(seq_len(ng), each = 3L))
x <- rnorm(ng * 3L)
u <- rnorm(ng, 0, 2)
d2 <- data.frame(y = rbinom(ng * 3L, 1, plogis(-0.5 + 0.5 * x + u[g])),
                 x = x, g = g, ntr = 1L)
d3 <- impsamp_data("poisson", "slope", "40x5", 2)
cases <- list(
  list(lab = "p124 binomial", form = bf(y | trials(ntr) ~ x + (1 | g)),
       fam = binomial(), d = d1, famname = "binomial", shape = "int"),
  list(lab = "hard bernoulli", form = bf(y ~ x + (1 | g)), fam = bernoulli(),
       d = d2, famname = "bernoulli", shape = "int"),
  list(lab = "poisson slope 40x5", form = bf(y ~ x + (1 + x | g)),
       fam = poisson(), d = d3, famname = "poisson", shape = "slope"))
pair_mcse <- function(a) {
  nd <- ncol(a)
  w <- exp(a - apply(a, 1L, max))
  half <- nd %/% 2L
  pm <- 0.5 * (w[, seq_len(half), drop = FALSE] +
                 w[, half + seq_len(half), drop = FALSE])
  m <- rowMeans(pm)
  v <- apply(pm, 1L, stats::var) / half
  sqrt(sum(v / m^2))
}
for (cc in cases) {
  for (N in c(64L, 200L, 512L)) {
    rows <- list()
    for (s in 1:40) {
      fit <- tryCatch(suppressWarnings(suppressMessages(frm(
        cc$form, family = cc$fam, data = cc$d, importance = N,
        control = frmtmb_control(importance_seed = s)))),
        error = function(e) e)
      if (inherits(fit, "error")) {
        rows[[s]] <- c(ll = NA, ex = NA, iid = NA, pair = NA)
        next
      }
      prop <- ns$imp_frozen_proposal(fit)
      a <- as.matrix(prop$io$amat(ns$imp_par_list(
        prop$template, fit$opt$par, prop$map)))
      rows[[s]] <- c(ll = as.numeric(logLik(fit)),
                     ex = impsamp_exact_at(cc$famname, cc$shape, cc$d, fit),
                     iid = fit$importance$mcse, pair = pair_mcse(a))
    }
    r <- as.data.frame(do.call(rbind, rows))
    e <- r$ll - r$ex
    cat(sprintf(paste0("%s N=%d fits %d: mean err %+.4f sd err %.4f | median ",
                       "iid MCSE %.4f pair MCSE %.4f | sd/iid %.2f sd/pair %.2f",
                       " | cover95 iid %.2f pair %.2f\n"),
                cc$lab, N, sum(!is.na(e)), mean(e, na.rm = TRUE),
                sd(e, na.rm = TRUE), median(r$iid, na.rm = TRUE),
                median(r$pair, na.rm = TRUE),
                sd(e, na.rm = TRUE) / median(r$iid, na.rm = TRUE),
                sd(e, na.rm = TRUE) / median(r$pair, na.rm = TRUE),
                mean(abs(e) < 1.96 * r$iid, na.rm = TRUE),
                mean(abs(e) < 1.96 * r$pair, na.rm = TRUE)))
  }
}
