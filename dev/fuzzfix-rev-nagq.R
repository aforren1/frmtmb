# Reviewer: what lme4, glmmTMB and GLMMadaptive do with adaptive
# quadrature over nested or several random-effect terms.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages({library(lme4); library(glmmTMB); library(GLMMadaptive)})
cat("lme4", format(packageVersion("lme4")), "glmmTMB",
    format(packageVersion("glmmTMB")), "GLMMadaptive",
    format(packageVersion("GLMMadaptive")), "\n")
set.seed(31)
ng <- 12; m <- 4
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
d$x <- rnorm(nrow(d))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$y <- rbinom(nrow(d), 1, plogis(0.3 + 0.5 * d$x + rnorm(ng, 0, 0.8)[d$ga] +
                                   rnorm(nlevels(ab), 0, 0.6)[ab]))
d$ab <- ab
show <- function(lbl, expr) {
  r <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  cat(sprintf("%-40s %s\n", lbl,
              if (inherits(r, "condition")) paste(class(r)[1], conditionMessage(r))
              else sprintf("ok logLik %.8f", as.numeric(logLik(r)))))
}
show("glmer (1|ga), nAGQ = 25",
     glmer(y ~ x + (1 | ga), d, binomial, nAGQ = 25))
show("glmer (1|ga/gb), nAGQ = 1",
     glmer(y ~ x + (1 | ga/gb), d, binomial, nAGQ = 1))
show("glmer (1|ga/gb), nAGQ = 5",
     glmer(y ~ x + (1 | ga/gb), d, binomial, nAGQ = 5))
show("glmer (1|ga) + (1|ab), nAGQ = 5",
     glmer(y ~ x + (1 | ga) + (1 | ab), d, binomial, nAGQ = 5))
show("glmer (x||ga), nAGQ = 5",
     glmer(y ~ x + (x || ga), d, binomial, nAGQ = 5))
cat("glmmTMB has no nAGQ argument:",
    !"nAGQ" %in% names(formals(glmmTMB::glmmTMB)), "\n")
show("glmmTMB (1|ga/gb) (Laplace)",
     glmmTMB(y ~ x + (1 | ga/gb), d, family = binomial))
show("GLMMadaptive random = ~ 1 | ga",
     mixed_model(y ~ x, random = ~ 1 | ga, data = d, family = binomial(),
                 nAGQ = 15))
show("GLMMadaptive random = ~ 1 | ga/gb",
     mixed_model(y ~ x, random = ~ 1 | ga/gb, data = d, family = binomial(),
                 nAGQ = 15))
show("GLMMadaptive random = list(~1|ga, ~1|ab)",
     mixed_model(y ~ x, random = list(~ 1 | ga, ~ 1 | ab), data = d,
                 family = binomial(), nAGQ = 15))
