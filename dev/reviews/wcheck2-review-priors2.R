# Review of lane wcheck2: do the prior errors that the 4985 test expects
# survive lane priors2's resolver? Runs the test's three priors on the
# unweighted model (no refusal) under priors2's library, and calls
# priors2's resolve_prior_input() the way the lane's new step does.
#   Rscript dev/reviews/wcheck2-review-priors2.R <lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")
set.seed(20261111)
d <- data.frame(obs = factor(seq_len(80)), x = stats::rnorm(80))
d$y <- 0.5 + 0.4 * d$x + stats::rnorm(80)
d$w <- 1.1
fx <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + x)
fu <- bf(y ~ x + (1 | obs), sigma ~ 1 + x)
P <- function(...) set_prior(..., dpar = "sigma")
pri <- list(
  coef_intercept = P("normal(0, 1)", class = "b", coef = "Intercept"),
  coef_z = P("normal(0, 1)", class = "b", coef = "z"),
  t_intercept = P("student_t(3, 0, 2.5)", class = "Intercept"),
  lb_intercept = P("normal(0, 1)", class = "Intercept", lb = -5))
fr <- q(frm(fx, data = d, dry_run = "frame"))
for (nm in names(pri)) {
  e1 <- tryCatch({ q(frm(fu, data = d, prior = pri[[nm]])); "fit" },
                 error = function(e) conditionMessage(e))
  e2 <- tryCatch({ ns$resolve_prior_input(list(frame = fr, spec = fr$spec),
                                          pri[[nm]]); "resolves" },
                 error = function(e) conditionMessage(e))
  e3 <- tryCatch({ q(frm(fx, data = d, prior = pri[[nm]])); "fit" },
                 error = function(e) conditionMessage(e))
  cat(sprintf("%-15s unweighted: %s\n%-15s resolver on the weighted frame: %s\n%-15s weighted frm(): %s\n",
              nm, substr(gsub("\\s+", " ", e1), 1, 100), "",
              substr(gsub("\\s+", " ", e2), 1, 100), "",
              substr(gsub("\\s+", " ", e3), 1, 100)))
}
cat("DONE\n")
