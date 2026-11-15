# Review of lane nlstop2, attack 3: what the huber() quadrature route
# changes downstream of the fit, and the model shapes it reaches.
# REV_ARM=before|after Rscript dev/reviews/nlstop2-review-huber.R
source("dev/reviews/nlstop2-review-common.R")
hd <- function(seed, n = 100, ng = 20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), length.out = n)),
                  x = rnorm(n))
  b <- rnorm(ng)[d$g]
  d$y <- 1 + 0.5 * d$x + 0.5 * b + rt(n, 3)
  d$h <- factor(rep(1:10, length.out = n))
  d$a <- factor((as.integer(d$g) - 1) %/% 4)  # g nested in a
  d$z <- rnorm(n)
  d$w <- runif(n, 0.5, 2)
  d$y2 <- rnorm(n) + 0.3 * b
  d$f <- factor(rep(c("p", "q"), length.out = n))
  d
}
d <- hd(4)
t0 <- proc.time()[[3]]
r <- tt(frm(y ~ x + (1 | g), family = huber(), data = d))
show("fit (1|g)", r)
fit <- r$v
cat("quadrature", isTRUE(fit$quadrature), "ll", llv(fit), "t",
    proc.time()[[3]] - t0, "\n")

# downstream methods on the default fit
show("summary", tt(summary(fit)))
show("confint wald", tt(confint(fit)))
show("vcov", tt(vcov(fit)))
show("fitted se", tt(fitted(fit, se.fit = TRUE)))
show("predict se", tt(predict(fit, se.fit = TRUE)))
show("predict newdata se", tt(predict(fit, newdata = d[1:5, ],
                                      se.fit = TRUE)))
show("predict brms-like", tt(predict(fit, ndraws = 50)))
show("posterior_epred", tt(brms::posterior_epred(fit, ndraws = 50)))
show("residuals osa", tt(residuals(fit, type = "osa")))
show("residuals response", tt(residuals(fit)))
show("dharma", tt(dharma_residuals(fit)))
show("ranef condVar", tt(ranef(fit, condVar = TRUE)))
show("VarCorr", tt(VarCorr(fit)))
show("vcov_cluster", tt(vcov_cluster(fit)))
show("cluster_scores", tt(cluster_scores(fit)))
show("simulate", tt(simulate(fit, nsim = 2)))
show("update", tt(update(fit, . ~ . + z)))
f2 <- tt(update(fit, . ~ . + z))$v
show("anova", tt(anova(fit, f2)))
show("AIC", tt(AIC(fit)))
show("loo", tt(loo(fit)))
show("influence", tt(influence(fit)))
show("drop1", tt(drop1(fit)))
show("hypothesis", tt(hypothesis(fit, "x = 0")))
show("conditional_effects", tt(conditional_effects(fit)))
show("confint profile", tt(confint(fit, method = "profile", parm = "x")))
show("frm_bootstrap", tt(frm_bootstrap(fit, nsim = 3)))
show("refit", tt(refit(fit, newresp = d$y + 0.1)))
show("allfit", tt(frm_allfit(fit)))

# model shapes the route reaches
shape <- function(tag, expr) {
  t0 <- proc.time()[[3]]
  r <- tt(expr)
  el <- proc.time()[[3]] - t0
  show(tag, r)
  if (!inherits(r$v, "tterr")) {
    cat("   quad", isTRUE(r$v$quadrature), "ll", llv(r$v), "code",
        r$v$opt$convergence, "t", round(el, 2), "\n")
  }
  invisible(r$v)
}
shape("nested (1|a/g)", frm(y ~ x + (1 | a / g), family = huber(),
                             data = d))
shape("nested (1|a)+(1|a:g)", frm(y ~ x + (1 | a) + (1 | a:g),
                                   family = huber(), data = d))
shape("crossed (1|g)+(1|h)", frm(y ~ x + (1 | g) + (1 | h),
                                  family = huber(), data = d))
shape("slope (0+x|g)", frm(y ~ x + (0 + x | g), family = huber(),
                            data = d))
shape("corr (1+x|g)", frm(y ~ x + (1 + x | g), family = huber(),
                           data = d))
shape("uncorr (1+x||g)", frm(y ~ x + (1 + x || g), family = huber(),
                              data = d))
shape("weights", frm(y | weights(w) ~ x + (1 | g), family = huber(),
                     data = d))
shape("sigma (1|g)", frm(bf(y ~ x, sigma ~ (1 | g)), family = huber(),
                         data = d))
shape("mu+sigma (1|g)", frm(bf(y ~ x + (1 | g), sigma ~ (1 | g)),
                            family = huber(), data = d))
shape("s(x)+(1|g)", frm(y ~ s(z) + (1 | g), family = huber(), data = d))
shape("by gr", frm(y ~ x + (1 | gr(g, by = f)), family = huber(),
                   data = d))
shape("mv huber+gaussian", frm(mvbf(bf(y ~ x + (1 | g)),
                                    bf(y2 ~ x + (1 | g)), rescor = FALSE),
                               family = list(huber(), gaussian()),
                               data = d))
shape("mv huber+huber shared", frm(mvbf(bf(y ~ x + (1 | p | g)),
                                        bf(y2 ~ x + (1 | p | g)),
                                        rescor = FALSE),
                                   family = list(huber(), huber()),
                                   data = d))
shape("mixture huber", frm(y ~ x + (1 | g),
                           family = mixture(huber(), huber()), data = d))
shape("student re", frm(y ~ x + (1 | gr(g, dist = "student")),
                        family = huber(), data = d))
shape("nl body", frm(bf(y ~ a0 + a1 * x, a0 ~ 1 + (1 | g), a1 ~ 1,
                        nl = TRUE), family = huber(), data = d))
shape("explicit NULL", frm(y ~ x + (1 | g), family = huber(), data = d,
                           quadrature = NULL))
shape("start given", frm(y ~ x + (1 | g), family = huber(), data = d,
                         start = list(beta = c(1, 0.5))))
shape("prior", frm(y ~ x + (1 | g), family = huber(), data = d,
                   prior = prior(normal(0, 1), class = "b")))
cat("DONE\n")
