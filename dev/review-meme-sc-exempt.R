# Review of lane meme, short check: the exemptions of the coef-bound
# refusal. Rscript dev/review-meme-sc-exempt.R > dev/review-meme-out/sc/exempt.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
q <- function(expr) suppressMessages(suppressWarnings(expr))
try_ <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1) cat("  ", r, "\n") else print(r)
  invisible(r)
}
fn <- bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)
set.seed(2)
dn <- data.frame(x = runif(60, 0, 3))
dn$y <- 2 * exp(-0.7 * dn$x) + rnorm(60, 0, 0.1)
run <- function(pr) {
  ds <- q(frm_sample(fn, data = dn, family = gaussian(),
                     prior = pr + set_prior("exponential(1)", class = "sigma"),
                     sample_prior = "yes", chains = 2, iter = 600, seed = 3,
                     refresh = 0))
  pd <- prior_draws(ds)
  m <- as.matrix(ds)
  list(prior_names = names(pd),
       prior_min = vapply(pd, min, 0),
       post_min = apply(m[, grep("^b_[ab]_", colnames(m)), drop = FALSE], 2,
                        min))
}
cat("#### nlpar intercepts, a ~ 1 and b ~ 1\n")
try_("frmtmb spelling: class Intercept, nlpar, lb (exempt)",
     run(set_prior("normal(1, 1)", class = "Intercept", nlpar = "a", lb = 0) +
           set_prior("normal(1, 1)", class = "Intercept", nlpar = "b",
                     lb = 0.4)))
try_("brms spelling, class row: class b, nlpar, lb (no coef)",
     run(set_prior("normal(1, 1)", class = "b", nlpar = "a", lb = 0) +
           set_prior("normal(1, 1)", class = "b", nlpar = "b", lb = 0.4)))
cat("\n#### theta on a toep() block\n")
set.seed(4)
dd <- data.frame(g = factor(rep(1:30, each = 4)), tim = factor(rep(1:4, 30)))
dd$y <- rnorm(120) + rnorm(30)[dd$g]
ft <- bf(y ~ 1 + toep(tim + 0 | g))
dp <- as.data.frame(default_prior(ft, data = dd))
print(dp[dp$class %in% c("theta", "sd", "cor"),
         c("prior", "class", "coef", "group", "lb", "ub")])
th <- dp$coef[dp$class == "theta" & nzchar(dp$coef)][1]
if (!is.na(th) && nzchar(th)) {
  try_(paste("frm_sample theta coef", th, "with a bound"), {
    ds <- frm_sample(ft, data = dd, family = gaussian(),
                     prior = set_prior("", class = "theta", coef = th,
                                       lb = -0.5, ub = 0.5),
                     chains = 1, iter = 200, seed = 1, refresh = 0)
    "sampled"
  })
}
cat("DONE\n")
