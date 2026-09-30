# Reviewer, lane sampfix, script 07: is a FULL draws object ever taken for
# laplace draws? draws_index() now reads every draw through
# draws_is_laplace(), so a false positive would read full draws in the
# laplace layout. Models with few random effects next to a constant
# dpar (mapped betad entries absent from the draws), a mixture (a
# derived weight column) and me() latent values. Real sampling, seed 3.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(12)
dd <- data.frame(g = factor(rep(1:2, each = 20)), x = rnorm(40))
dd$y <- 0.5 + 0.4 * dd$x + c(-0.3, 0.3)[dd$g] + rnorm(40, 0, 0.5)
dd$ym <- ifelse(runif(40) < 0.5, rnorm(40, -2), rnorm(40, 2)) + c(-0.3, 0.3)[dd$g]
dd$xe <- dd$x + rnorm(40, 0, 0.2); dd$sdx <- 0.2
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
smp <- function(fit, ...) suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3, ...)))
cases <- list(
  sigma_const = function() q(bf(y ~ x + (1 | g), sigma = 0.5), family = gaussian(), data = dd),
  mix_2lev = function() q(bf(ym ~ 1 + (1 | g)), family = mixture(gaussian(), gaussian()), data = dd),
  mix_theta_const = function() q(bf(ym ~ 1 + (1 | g), theta1 = 0), family = mixture(gaussian(), gaussian()), data = dd),
  me_term = function() q(bf(y ~ me(xe, sdx)), family = gaussian(), data = dd)
)
for (nm in names(cases)) {
  fit <- tryCatch(cases[[nm]](), error = function(e) e)
  if (inherits(fit, "error")) { cat(nm, "FIT ERROR", conditionMessage(fit), "\n"); next }
  for (lap in c(FALSE, TRUE)) {
    ds <- tryCatch(smp(fit, laplace = lap), error = function(e) e)
    if (inherits(ds, "error")) { cat(sprintf("%-16s laplace=%s SAMPLE ERROR %s\n", nm, lap, substr(conditionMessage(ds), 1, 100))); next }
    il <- frmtmb.sample:::draws_is_laplace(ds)
    ep <- tryCatch(posterior_epred(ds, ndraws = 5, re_formula = NA), error = function(e) e)
    cat(sprintf("%-16s laplace=%-5s ncol %d  draws_is_laplace %-5s  correct %s  epred(re NA) %s\n",
                nm, lap, ncol(ds$draws), il, identical(il, lap && nm != "none"),
                if (inherits(ep, "error")) paste("ERROR", substr(conditionMessage(ep), 1, 60)) else
                  sprintf("finite %s", all(is.finite(ep)))))
  }
}
cat("DONE\n")
