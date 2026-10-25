# Reviewer: a third-party option term written from the help page of
# frmtmb_register_aterm() alone.
#   Rscript dev/review-bhaz-toy.R > dev/review-bhaz-log/toy.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
# A normal family with a known sd that the term sets: y | sdk(2) ~ x.
knownsd <- function(sd = 1) {
  force(sd)
  fam <- frmtmb_family(
    "knownsd", dpars = "mu", links = list(mu = "identity"),
    accepts_aterms = c("sdk", "sdm"),
    lpdf = function(y, dpars, aterms, extra) {
      -0.5 * ((y - dpars$mu) / sd)^2 - log(sd) - 0.5 * log(2 * pi)
    })
  fam$knownsd_sd <- sd
  fam
}
set.seed(4); dk <- data.frame(x = rnorm(100)); dk$y <- 1 + dk$x +
  rnorm(100, 0, 2)
ref <- function(s) sum(stats::dnorm(dk$y, fitted(lm(y ~ x, dk)), s,
                                    log = TRUE))
f0 <- tryCatch(frm(bf(y ~ x), family = knownsd(2), data = dk),
               error = function(e) e)
cat("family alone, sd = 2:", if (inherits(f0, "error"))
  conditionMessage(f0) else sprintf("logLik %.8f ref %.8f",
                                    as.numeric(logLik(f0)), ref(2)), "\n")
frmtmb_local_registries()
calls <- 0L
frmtmb_register_aterm("sdk", option = function(family, call, env) {
  calls <<- calls + 1L
  knownsd(eval(call[[2L]], env))
})
ft <- tryCatch(frm(bf(y | sdk(2) ~ x), family = knownsd(), data = dk),
               error = function(e) e)
if (inherits(ft, "error")) {
  cat("option term: ERROR", conditionMessage(ft), "\n")
} else {
  cat(sprintf("option term sdk(2): logLik %.8f ref %.8f; calls %d\n",
              as.numeric(logLik(ft)), ref(2), calls))
  calls <- 0L
  u <- update(ft, . ~ . + I(x^2))
  cat("update(+ I(x^2)) calls", calls, "\n")
  calls <- 0L
  lp <- frm_linpred(ft, newdata = dk[1:3, ])
  cat("frm_linpred(newdata) calls", calls, "\n")
  gp <- get_prior(bf(y | sdk(2) ~ x), family = knownsd(), data = dk)
  cat("get_prior calls", calls, "\n")
}
# The help says option() must give the same family for the same family
# and call, and that update() parses again from the family as written.
# A multiplying option shows which family update() hands it.
frmtmb_register_aterm("sdm", option = function(family, call, env) {
  knownsd(family$knownsd_sd * eval(call[[2L]], env))
})
fm <- frm(bf(y | sdm(2) ~ x), family = knownsd(1), data = dk)
u <- update(fm, . ~ . + I(x^2))
u2 <- update(fm, newdata = dk[1:80, ])
cat(sprintf("non-idempotent sdm(2) from knownsd(1): fit %.8f (ref sd 2 %.8f)\n",
            as.numeric(logLik(fm)), ref(2)))
cat("family(fit)$knownsd_sd", family(fm)$knownsd_sd,
    " family(update)$knownsd_sd", family(u)$knownsd_sd,
    " family(update newdata)$knownsd_sd", family(u2)$knownsd_sd, "\n")
for (nm in c("sdk", "sdk()")) {
  r <- tryCatch(frm_compat(nm), error = function(e) conditionMessage(e))
  cat("frm_compat(", nm, "):\n"); print(r)
}
# a term on a family that declares no allow-list
r <- tryCatch(frm(bf(y | sdk(2) ~ x), family = gaussian(), data = dk),
              error = function(e) conditionMessage(e))
cat("sdk() on gaussian():", if (is.character(r)) r else "FIT", "\n")
# the doc example: bw(5) sets family$width
frmtmb_register_aterm("bw", option = function(family, call, env) {
  family$width <- eval(call[[2L]], env)
  family
})
r <- tryCatch(frm(bf(y | bw(5) ~ x), family = gaussian(), data = dk),
              error = function(e) conditionMessage(e))
cat("doc example bw(5) on gaussian():", if (is.character(r)) r else
  paste("FIT; family$width", family(r)$width), "\n")
