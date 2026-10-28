# Follow-ups to dev/priordraws-evidence.R, on the lane's build.
#   Part 5: brms's own prior draws of model B against the exact prior
#           (one-sample KS), to read part 2's two-sample p-values.
#   Part 6: sample_prior = "only" on brms and frmtmb at two more seeds,
#           each side's quantiles against the exact prior quantiles,
#           z-scores with posterior::mcse_quantile().
# Usage: Rscript dev/priordraws-evidence2.R [part ...]
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
loadNamespace("brms")
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
parts <- commandArgs(TRUE)
if (!length(parts)) parts <- c("5", "6")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
brmq <- function(...) quiet(brms::brm(..., refresh = 0, cores = 1))
half_t <- function(q) pmax(2 * pt(q / 2.5, 3) - 1, 0)

if ("5" %in% parts) {
  cat("\n#### PART 5: brms's own prior draws against the exact prior\n")
  pB <- c(brms::prior(normal(0, 1), class = b),
          brms::prior(normal(0, 2), class = b, coef = z),
          brms::prior(normal(1, 3), class = Intercept),
          brms::prior(exponential(2), class = sigma),
          brms::prior(lkj(2), class = cor),
          brms::prior(cauchy(0, 1), class = sd, group = g, coef = x))
  for (s in c(11, 12)) {
    fb <- brmq(y ~ x + z + (1 + x | g), data = d, prior = pB,
               sample_prior = "yes", chains = 2, iter = 4000, seed = s)
    pb <- brms::prior_draws(fb)
    ref <- list(
      Intercept = function(q) pnorm(q, 1, 3), b_x = pnorm,
      b_z = function(q) pnorm(q, 0, 2), sigma = function(q) pexp(q, 2),
      sd_g__Intercept = half_t,
      sd_g__x = function(q) pmax(2 * pcauchy(q) - 1, 0),
      cor_g = function(q) pbeta((q + 1) / 2, 2, 2))
    for (cn in names(ref)) {
      ks <- ks.test(pb[[cn]], ref[[cn]])
      cat(sprintf("  seed %d brms %-16s n=%d KS D=%.4f p=%.3f\n", s, cn,
                  length(pb[[cn]]), ks$statistic, ks$p.value))
    }
  }
}

if ("6" %in% parts) {
  cat("\n#### PART 6: sample_prior = \"only\", each side against the exact prior\n")
  exact <- list(
    b_x = function(p) qnorm(p),
    sigma = function(p) 2.5 * qt((1 + p) / 2, 3),
    sd_g__Intercept = function(p) 2.5 * qt((1 + p) / 2, 3))
  for (s in c(12, 13)) {
    fb <- brmq(y ~ x + (1 | g), data = d,
               prior = brms::prior(normal(0, 1), class = b),
               sample_prior = "only", chains = 4, iter = 4000, seed = s)
    ds <- quiet(frm_sample(frmtmb::bf(y ~ x + (1 | g)), data = d,
                           family = gaussian(),
                           prior = frmtmb::set_prior("normal(0, 1)",
                                                     class = "b"),
                           sample_prior = "only", chains = 4, iter = 4000,
                           seed = s, refresh = 0))
    for (side in c("brms", "frmtmb")) {
      m <- if (side == "brms") as.matrix(fb) else as.matrix(ds)
      for (v in names(exact)) {
        a <- matrix(m[, v], ncol = 4)
        zs <- vapply(c(0.05, 0.5, 0.95), function(p) {
          (quantile(a, p) - exact[[v]](p)) /
            posterior::mcse_quantile(a, probs = p)
        }, 1)
        cat(sprintf("  seed %d %-6s %-16s q05 %7.3f q50 %7.3f q95 %7.3f | exact %7.3f %7.3f %7.3f | z %6.2f %6.2f %6.2f | ess_bulk %6.0f\n",
                    s, side, v, quantile(a, .05), quantile(a, .5),
                    quantile(a, .95), exact[[v]](.05), exact[[v]](.5),
                    exact[[v]](.95), zs[1], zs[2], zs[3],
                    posterior::ess_bulk(a)))
      }
    }
  }
}
