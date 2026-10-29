# Release fix R1: how frm() (core) resolves the pair of spellings of a
# ~ 1's intercept, set_prior(class = "Intercept", nlpar = "a") (pA) and
# brms's class row set_prior(class = "b", nlpar = "a") (pB), with and
# without bounds, in both orders: the resolved density and bounds on
# a_(Intercept) and frm()'s MAP. frm_sample() is made to resolve the pair
# the same way (dev/rel078-release-fix.R). Data as the reviewer's
# dev/review-meme-sc2.R (seed 2).
#   Rscript dev/rel078-r1-frm.R [lib] > dev/rel078-log/r1-frm.txt
lib <- commandArgs(TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r16"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("lib:", lib, "\n")
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressMessages(suppressWarnings(e))
set.seed(2)
d <- data.frame(x = runif(60, 0, 3), w = rnorm(60))
d$y <- 2 * exp(-0.7 * d$x) + rnorm(60, 0, 0.1)
f1 <- bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)
ents <- function(pr) {
  fit <- q(frm(f1, data = d, dry_run = "objective", prior = pr))
  ri <- q(frmtmb:::resolve_prior_input(fit, pr))
  pt <- fit$frame$par_template
  e <- vapply(ri$entries, function(e) {
    nm <- frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
    paste0(paste(nm, collapse = ","), ":", e$dist$kind, "(",
           paste(unlist(e$dist[-1]), collapse = ","), ")")
  }, "")
  cat(" entries:", e, "\n lower:", format(unlist(ri$lower)), "\n upper:",
      format(unlist(ri$upper)), "\n")
}
est <- function(pr) {
  r <- tryCatch({
    fit <- q(frm(f1, data = d, prior = pr))
    sprintf("%.4f", fit$estimates$beta)
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(" frm beta:", r, "\n")
}
cases <- list(
  "A dens + B dens" = c("normal(5, 0.1)", NA, "normal(0, 0.1)", NA),
  "A dens + B dens lb 6" = c("normal(5, 0.1)", NA, "normal(0, 0.1)", 6),
  "A dens lb 6 + B dens" = c("normal(5, 0.1)", 6, "normal(0, 0.1)", NA),
  "A dens lb 1 + B dens lb 6" = c("normal(5, 0.1)", 1, "normal(0, 0.1)", 6),
  "A dens lb 6 + B dens lb 1" = c("normal(5, 0.1)", 6, "normal(0, 0.1)", 1),
  "A empty lb 6 + B dens" = c("", 6, "normal(0, 0.1)", NA),
  "A dens + B empty lb 6" = c("normal(5, 0.1)", NA, "", 6),
  "A empty + B dens" = c("", NA, "normal(0, 0.1)", NA)
)
for (nm in names(cases)) {
  cc <- cases[[nm]]
  pA <- set_prior(cc[1], class = "Intercept", nlpar = "a",
                  lb = as.numeric(cc[2]))
  pB <- set_prior(cc[3], class = "b", nlpar = "a", lb = as.numeric(cc[4]))
  for (ord in c("AB", "BA")) {
    pr <- if (ord == "AB") pA + pB else pB + pA
    cat("==", nm, ord, "\n")
    tryCatch(ents(pr), error = function(e) cat(" ents ERROR",
                                                conditionMessage(e), "\n"))
    est(pr)
  }
}
