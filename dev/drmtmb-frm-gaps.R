# Models drmTMB fits that frmtmb refuses, each run in both packages.
# The first three families are brms families, so these are brms-parity
# gaps in frmtmb; the tree input is a drmTMB convenience brms lacks too.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n))
d$ypos <- rnbinom(n, mu = exp(1 + 0.3 * d$x), size = 2) + 1L
d$yh <- ifelse(runif(n) < 0.3, 0L, d$ypos)
r <- runif(n)
d$p <- ifelse(r < 0.1, 0, ifelse(r < 0.2, 1, rbeta(n, 4, 4)))
tree <- ape::rcoal(10)
d$sp <- factor(sample(tree$tip.label, n, TRUE))
d$yg <- rnorm(n)
run <- function(cl) {
  cat("\n", paste(deparse(cl, width.cutoff = 70), collapse = "\n"), "\n")
  r <- tryCatch(suppressWarnings(eval(cl)), error = function(e) e)
  msg <- if (inherits(r, "error")) {
    paste("REFUSED:", substr(conditionMessage(r), 1, 300))
  } else {
    paste("FITTED: logLik", format(as.numeric(logLik(r)), digits = 10),
          " convergence", r$opt$convergence)
  }
  cat(msg, "\n")
}
dbf <- drmTMB::drm_formula
run(quote(drmTMB::drmTMB(dbf(yh ~ x, sigma ~ 1, hu ~ 1),
                         family = drmTMB::truncated_nbinom2(), data = d)))
run(quote(frm(bf(yh ~ x, hu ~ 1), family = "hurdle_negbinomial",
              data = d)))
run(quote(drmTMB::drmTMB(dbf(ypos ~ x, sigma ~ 1),
                         family = drmTMB::truncated_nbinom2(), data = d)))
run(quote(frm(bf(ypos | trunc(lb = 1) ~ x), family = negbinomial(),
              data = d)))
run(quote(drmTMB::drmTMB(dbf(p ~ x, sigma ~ 1, zoi ~ 1, coi ~ 1),
                         family = drmTMB::zero_one_beta(), data = d)))
run(quote(frm(bf(p ~ x, zoi ~ 1, coi ~ 1),
              family = "zero_one_inflated_beta", data = d)))
run(quote(drmTMB::drmTMB(dbf(yg ~ x + phylo(1 | sp, tree = tree),
                             sigma ~ 1), family = gaussian(), data = d)))
run(quote(frm(bf(yg ~ x + (1 | gr(sp, cov = tree))), family = gaussian(),
              data = d, data2 = list(tree = tree))))
