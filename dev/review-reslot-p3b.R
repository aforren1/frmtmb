# Reviewer probe 3b: default_prior() on an rr() model, both arms. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
set.seed(1)
G <- 8L; n <- 12L
dd <- data.frame(g = factor(rep(seq_len(G), each = n)), t = rep(seq_len(n), G))
dd$x <- rnorm(nrow(dd)); dd$tf <- factor((dd$t - 1) %% 3)
dd$y <- 1 + 0.5 * dd$x + rnorm(G, 0, 0.7)[dd$g] + rnorm(nrow(dd))
for (route in c("fit", "sample")) {
  r <- tryCatch(default_prior(bf(y ~ x + rr(0 + tf | g, d = 1)), data = dd,
                              route = route), error = function(e) e)
  cat(route, ":", if (inherits(r, "error")) conditionMessage(r) else
    paste(r$class, r$coef, collapse = "; "), "\n")
}
r <- tryCatch(validate_prior(set_prior("normal(0, 1)", class = "sd"),
                             bf(y ~ x + rr(0 + tf | g, d = 1)), data = dd),
              error = function(e) e)
cat("validate sd:", if (inherits(r, "error")) conditionMessage(r) else "ok", "\n")
fit <- tryCatch(frm(bf(y ~ x + rr(0 + tf | g, d = 1)), data = dd,
                    prior = set_prior("normal(0, 1)", class = "b")),
                error = function(e) e)
cat("frm with prior:", if (inherits(fit, "error")) conditionMessage(fit) else "ok", "\n")
