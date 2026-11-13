# Review of lane boundfam: do the new single-seed assertions hold over
# 20 ulp perturbations of the data (lane-env "Tests")? Lane build only.
.libPaths(c("C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/reviews/boundfam-review-fp.R")
eval(parse(text = src[grep("^conds <- ", src):(grep("^ok <- ", src) - 1L)]))
tally <- function(lab, mk) {
  out <- character(0)
  for (k in 0:19) {
    r <- conds(mk(k))
    lost <- frmtmb:::sdr_of(r$value)$se_lost
    key <- paste0("lost=", paste(names(lost), lost, sep = ":", collapse = ","),
                  " nW=", length(r$w), " nM=", length(r$m),
                  " Mclass=", paste(sub("^\\[([^]]*)\\].*", "\\1", r$m), collapse = ","))
    out <- c(out, key)
  }
  cat("\n## ", lab, "\n", sep = "")
  tb <- table(out)
  for (i in seq_along(tb)) cat(sprintf("  %2d x %s\n", tb[[i]], names(tb)[i]))
}
set.seed(2)
dnb <- data.frame(x = rnorm(200))
dnb$y <- rpois(200, exp(1 + 0.4 * dnb$x))
tally("negbinomial poisson data seed 2 (test-se-limit)", function(k) {
  d <- dnb; d$x <- d$x * (1 + 2e-15 * k)
  frm(y ~ x, family = negbinomial(), data = d)
})
for (s in c(1L, 3L)) {
  set.seed(s)
  dc <- data.frame(x = rnorm(60))
  dc$y <- 1 + dc$x + rnorm(60)
  dc$cen <- "right"
  tally(paste0("cens all right seed ", s, " (test-se-limit)"), function(k) {
    d <- dc; d$y <- d$y * (1 + 2e-15 * k)
    frm(y | cens(cen) ~ x, data = d)
  })
}
dm <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30, 20, 6,
                       44, 44, 27, 25, 17, 14, 27, 35, 30), n = 45L)
tally("Malingering_2 (test-se-limit), start moved by ulps", function(k) {
  frm(bf(k | trials(n) ~ 1), family = mixture(beta_binomial, beta_binomial),
      data = dm, start = list(beta = qlogis(c(0.55, 0.95)) * (1 + 2e-15 * k),
                              betad = c(log(c(15, 15)), 0)))
})
set.seed(41)
du <- data.frame(x = rnorm(120))
du$y <- 1 + 0.5 * du$x + (runif(120) - 0.5) * 3.4
tally("student uniform errors (test-diagnostics-ux)", function(k) {
  d <- du; d$y <- d$y * (1 + 2e-15 * k)
  frm(bf(y ~ x) + student(), data = d)
})
set.seed(77101)
dcy <- data.frame(x = rnorm(200))
dcy$y <- 1 + 0.5 * dcy$x + rcauchy(200)
tally("student cauchy errors (test-diagnostics-ux)", function(k) {
  d <- dcy; d$y <- d$y * (1 + 2e-15 * k)
  frm(bf(y ~ x) + student(), data = d)
})
