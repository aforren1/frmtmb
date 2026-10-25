## The old design-column spelling of a coef, and brms's spelling, in
## every place a user writes one: set_prior(), prior_string(),
## validate_prior(), prior_summary() of an frm() fit and of
## frm_sample() draws. A tight normal(3, 0.001) prior shows which
## coefficient the fit actually moved. Seed 31.
##   Rscript dev/review-brmsnames-prior.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(library(frmtmb.sample))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(31)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), xp = runif(n, 0.5, 2),
                k = factor(sample(1:3, n, replace = TRUE)))
d$y <- 1 + d$x + 0.3 * d$z^2 + 0.5 * log(d$xp) + rnorm(n)
form <- bf(y ~ poly(x, 2) + I(z^2) + log(xp) + factor(k))
fit0 <- q(frm(form, data = d))
cat("fixef rows:", rownames(fixef(fit0)), "\n")
pairs <- list(c("poly(x, 2)1", "polyx21"), c("I(z^2)", "IzE2"),
              c("log(xp)", "logxp"), c("factor(k)2", "factork2"))
row_of <- function(tab, val) {
  tab <- as.data.frame(tab)
  u <- tab[tab$source == "user", , drop = FALSE]
  paste0(nrow(tab), " rows; user rows: ",
         paste(sprintf("[%s|%s]", u$coef, u$prior), collapse = " "))
}
for (pr in pairs) {
  brmsname <- pr[2]
  for (spell in pr) {
    tag <- sprintf("%-12s via %-12s", brmsname, shQuote(spell))
    p <- tryCatch(set_prior("normal(3, 0.001)", class = "b", coef = spell),
                  error = function(e) e)
    if (inherits(p, "error")) {
      cat(tag, "set_prior ERROR", conditionMessage(p), "\n")
      next
    }
    ps <- tryCatch(prior_string("normal(3, 0.001)", class = "b",
                                coef = spell), error = function(e) e)
    v <- tryCatch(q(validate_prior(p, form, data = d)),
                  error = function(e) e)
    f <- tryCatch(q(frm(form, data = d, prior = p)), error = function(e) e)
    est <- if (inherits(f, "error")) conditionMessage(f) else {
      fe <- fixef(f)[, 1]
      sprintf("moved: %s", paste(names(fe)[abs(fe - 3) < 0.01],
                                 collapse = ","))
    }
    ps_sum <- if (inherits(f, "error")) "" else row_of(prior_summary(f))
    cat(tag, "\n   prior_string:", if (inherits(ps, "error"))
      conditionMessage(ps) else "ok",
      "\n   validate_prior:", if (inherits(v, "error"))
        conditionMessage(v) else row_of(v),
      "\n   frm fit:", est,
      "\n   prior_summary(fit):", ps_sum, "\n")
  }
}
# frm_sample(prior =) with the old spelling, on the fit route
p <- set_prior("normal(3, 0.001)", class = "b", coef = "poly(x, 2)1")
ds <- tryCatch(q(frm_sample(fit0, prior = p, chains = 1, iter = 60,
                            warmup = 30, refresh = 0, seed = 1)),
               error = function(e) e)
if (inherits(ds, "error")) {
  cat("frm_sample old spelling ERROR:", conditionMessage(ds), "\n")
} else {
  m <- colMeans(as.matrix(ds)[, grep("^b_", variables(ds)), drop = FALSE])
  cat("frm_sample old spelling, draws mean near 3:",
      names(m)[abs(m - 3) < 0.05], "\n")
  cat("prior_summary(draws):", row_of(prior_summary(ds)), "\n")
}
p <- set_prior("normal(3, 0.001)", class = "b", coef = "IzE2")
ds <- tryCatch(q(frm_sample(fit0, prior = p, chains = 1, iter = 60,
                            warmup = 30, refresh = 0, seed = 1)),
               error = function(e) e)
if (inherits(ds, "error")) {
  cat("frm_sample brms spelling ERROR:", conditionMessage(ds), "\n")
} else {
  m <- colMeans(as.matrix(ds)[, grep("^b_", variables(ds)), drop = FALSE])
  cat("frm_sample brms spelling, draws mean near 3:",
      names(m)[abs(m - 3) < 0.05], "\n")
  cat("prior_summary(draws):", row_of(prior_summary(ds)), "\n")
}
