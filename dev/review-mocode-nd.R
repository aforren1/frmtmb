## Fix round 1 check: new-data refusals of mo() (fractional, NA, factor,
## character, out of range), marginaleffects slopes and comparisons,
## emmeans, mo(mi()). Seed 5 data as in dev/review-mocode-frm.R.
##   Rscript dev/review-mocode-nd.R lane
.libPaths(c("C:/Users/adf44/source/r/wt-mocode-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb lib:", dirname(system.file(package = "frmtmb")), "\n")
run <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(tryCatch(suppressMessages(expr),
                                    error = function(e) e),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("%-30s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", gsub("\n", " ", substr(conditionMessage(r), 1, 230)))
    else paste(format(signif(as.numeric(unlist(r)), 7)), collapse = " ")))
  if (length(w)) cat(sprintf("%-30s   warning: %s\n", "",
                             substr(unique(w)[1L], 1, 120)))
}
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fz <- suppressWarnings(suppressMessages(frm(y ~ mo(k) + z, data = d)))
nd <- function(x) data.frame(k = x, z = 0)
for (x in list(0.6, 2.5, 4.4, NA_integer_, c(2L, NA), 0L, 5L, 1L, 4L)) {
  run(paste("fitted k =", paste(x, collapse = ",")),
      fitted(fz, newdata = nd(x))[, 1])
}
run("fitted k = ordered factor", fitted(fz, newdata = nd(
  factor(2, levels = 1:4, ordered = TRUE)))[, 1])
run("fitted k = factor", fitted(fz, newdata = nd(factor(2))))
run("fitted k = '2'", fitted(fz, newdata = nd("2")))
run("fitted k = TRUE", fitted(fz, newdata = nd(TRUE)))
run("predict k = 2.5", predict(fz, newdata = nd(2.5)))
run("simulate k = 2.5", simulate(fz, nsim = 1, newdata = nd(2.5)))
run("avg_slopes(k)", marginaleffects::avg_slopes(fz, variables = "k")$estimate)
run("avg_slopes(z)", marginaleffects::avg_slopes(fz, variables = "z")$estimate)
run("avg_comparisons(k = c(1, 4))", {
  r <- marginaleffects::avg_comparisons(fz, variables = list(k = c(1, 4)))
  c(r$estimate, r$std.error)
})
run("avg_comparisons(k) default", {
  r <- marginaleffects::avg_comparisons(fz, variables = "k")
  c(r$estimate, r$std.error)
})
run("emmeans(~ z)", summary(emmeans::emmeans(fz, ~ z))$emmean)
run("emmeans(~ z, at k = 3)",
    summary(emmeans::emmeans(fz, ~ z, at = list(k = 3)))$emmean)
run("emmeans(~ k, at k = 1:4)",
    summary(emmeans::emmeans(fz, ~ k, at = list(k = 1:4)))$emmean)
run("emmeans cov.reduce", summary(emmeans::emmeans(
  fz, ~ z, cov.reduce = function(v) round(median(v))))$emmean)
di <- d
di$k[c(3, 7)] <- NA
run("mo(mi(k))", frm(bf(y ~ mo(mi(k))) + bf(k | mi() ~ 1) +
                       set_rescor(FALSE), data = di))
run("mo(me(k, 0.1))", frm(y ~ mo(me(k, 0.1)), data = d))
run("ref_grid cov.reduce fun: z", emmeans::ref_grid(
  fz, cov.reduce = function(v) round(median(v)))@grid$z)
run("emmeans cov.reduce list(k)", summary(emmeans::emmeans(
  fz, ~ z, cov.reduce = list(k = function(v) round(median(v)))))$emmean)
run("ref_grid cov.reduce list: z", emmeans::ref_grid(
  fz, cov.reduce = list(k = function(v) round(median(v))))@grid$z)
r <- tryCatch(marginaleffects::avg_slopes(fz, variables = "k"),
              error = function(e) conditionMessage(e))
cat("avg_slopes(k) full message:\n", r, "\n")
