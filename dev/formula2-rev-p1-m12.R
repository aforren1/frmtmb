# Reviewer, punch round 1: M1 (hmm equations) and M2 (re_formula with
# cmc = FALSE) on shapes of the reviewer's own.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
sh <- function(x) {
  if (inherits(x, "ERR")) paste("ERR:", substr(x, 1, 140)) else "ok"
}
cat("== M1: three-state hmm, tr12 = \"tr13\" (seed 4002)\n")
set.seed(4002)
G <- 25; Tn <- 30
P <- rbind(c(.8, .1, .1), c(.1, .8, .1), c(.1, .1, .8))
z <- unlist(lapply(seq_len(G), function(g) {
  s <- integer(Tn); s[1] <- sample(1:3, 1)
  for (t in 2:Tn) s[t] <- sample(1:3, 1, prob = P[s[t - 1], ])
  s
}))
dd <- data.frame(id = rep(seq_len(G), each = Tn), t = rep(seq_len(Tn), G))
dd$y <- rnorm(nrow(dd), c(-2, 1, 4)[z], 0.7)
h0 <- q(frm(bf(y ~ 1), family = hmm(K = 3, gaussian(), time = t,
                                    group = id), data = dd))
h1 <- q(frm(bf(y ~ 1, tr12 = "tr13"),
            family = hmm(K = 3, gaussian(), time = t, group = id),
            data = dd))
cat("fits:", sh(h0), sh(h1), "\n")
if (!inherits(h1, "ERR")) {
  cat("df free", attr(logLik(h0), "df"), "equated", attr(logLik(h1), "df"),
      "logLik", as.numeric(logLik(h0)), as.numeric(logLik(h1)), "\n")
  e <- h1$estimates$betad
  v <- q(variables(h1)); cat("variables:", if (inherits(v, "ERR")) v else v,
                             "\n")
  s <- q(summary(h1)); cat("summary:", sh(s), "\n")
  if (!inherits(s, "ERR")) print(s)
  cat("hypothesis tr13_Intercept = 0:",
      sh(q(hypothesis(h1, "tr13_Intercept = 0"))), "\n")
  cat("hypothesis tr12_Intercept = 0:",
      sh(q(hypothesis(h1, "tr12_Intercept = 0"))), "\n")
  cat("hypothesis sigma1 = sigma2 (class NULL):",
      sh(q(hypothesis(h1, "sigma1 = sigma2", class = NULL))), "\n")
  cat("fixef:", rownames(fixef(h1)), "\n")
  cat("coef:", sh(q(coef(h1))), " frm_linpred tr12 == tr13:",
      isTRUE(all.equal(q(frm_linpred(h1, dpar = "tr12")),
                       q(frm_linpred(h1, dpar = "tr13")))), "\n")
}
cat("== M1 regression: mixture sigma1 = sigma2 still lists sigma1\n")
set.seed(11)
dm <- data.frame(x = rnorm(300)); k <- rbinom(300, 1, 0.4)
dm$y <- ifelse(k == 1, rnorm(300, 3), rnorm(300, -1))
fm <- frm(bf(y ~ x, sigma1 = "sigma2"), family = mixture(gaussian(),
                                                         gaussian()), data = dm)
cat("variables:", variables(fm), "\n")
cat("coef names:", names(coef(fm)), "\n")

cat("== M2: re_formula on cmc = FALSE fits (seed 51)\n")
set.seed(51)
n <- 360
d <- data.frame(f = factor(sample(c("a", "b", "c"), n, TRUE)),
                x = rnorm(n), h = factor(rep(1:18, each = 20)),
                h2 = factor(sample(1:9, n, TRUE)))
u <- matrix(rnorm(18 * 3, 0, 0.6), 18)
d$y <- rnorm(n, as.numeric(d$f) + d$x + (d$f == "b") * u[d$h, 1] +
               (d$f == "c") * u[d$h, 2] + d$x * u[d$h, 3] +
               rnorm(9, 0, 0.5)[d$h2], 0.6)
fit <- q(frm(bf(y ~ x + (0 + f + x | h) + (1 | h2), cmc = FALSE), data = d))
cat("fit:", sh(fit), "\n")
nd <- d[1:10, ]
fe <- fixef(fit)[, "Estimate"]
re <- ranef(fit)
rh <- re$h[as.character(nd$h), , drop = FALSE]
eta_fix <- fe["Intercept"] + fe["x"] * nd$x
eta_h <- rh[, "fb"] * (nd$f == "b") + rh[, "fc"] * (nd$f == "c") +
  rh[, "x"] * nd$x
p_h <- q(frm_linpred(fit, newdata = nd, re_formula = ~ (0 + f + x | h)))
cat("re_formula ~ (0 + f + x | h), max |diff| / sd(y) to hand:",
    if (inherits(p_h, "ERR")) p_h else
      max(abs(unname(p_h) - unname(eta_fix + eta_h))) / sd(d$y), "\n")
eta_f <- rh[, "fb"] * (nd$f == "b") + rh[, "fc"] * (nd$f == "c")
p_f <- q(frm_linpred(fit, newdata = nd, re_formula = ~ (0 + f | h)))
cat("re_formula ~ (0 + f | h) (a subset), max |diff| / sd(y):",
    if (inherits(p_f, "ERR")) p_f else
      max(abs(unname(p_f) - unname(eta_fix + eta_f))) / sd(d$y), "\n")
p_x <- q(frm_linpred(fit, newdata = nd, re_formula = ~ (0 + x | h)))
cat("re_formula ~ (0 + x | h), max |diff| / sd(y):",
    if (inherits(p_x, "ERR")) p_x else
      max(abs(unname(p_x) - unname(eta_fix + rh[, "x"] * nd$x))) / sd(d$y),
    "\n")
p_bad <- q(frm_linpred(fit, newdata = nd, re_formula = ~ (1 + f | h)))
cat("re_formula ~ (1 + f | h) (not in the fit):", sh(p_bad), "\n")
cat("-- mu with cmc = TRUE and sigma with cmc = FALSE, same term text\n")
fit2 <- q(frm(bf(y ~ x + (0 + f | h)) +
                lf(sigma ~ 1 + (0 + f | h), cmc = FALSE), data = d))
cat("fit2:", sh(fit2), "\n")
if (!inherits(fit2, "ERR")) {
  cn <- lapply(fit2$re_blocks, function(b) lapply(b$components, `[[`,
                                                  "cnms"))
  print(cn)
  a <- q(fitted(fit2, newdata = nd, dpar = "sigma",
                re_formula = ~ (0 + f | h)))
  b <- q(fitted(fit2, newdata = nd, dpar = "sigma"))
  cat("sigma re_formula = its term == full:",
      if (inherits(a, "ERR")) a else isTRUE(all.equal(a, b)), "\n")
  a2 <- q(fitted(fit2, newdata = nd, dpar = "mu",
                 re_formula = ~ (0 + f | h)))
  b2 <- q(fitted(fit2, newdata = nd, dpar = "mu"))
  cat("mu re_formula = its term == full:",
      if (inherits(a2, "ERR")) a2 else isTRUE(all.equal(a2, b2)), "\n")
}
cat("== cmc = FALSE with a level-read structure over f:x1 (no dropped level)\n")
set.seed(41)
dc <- data.frame(g = factor(rep(1:30, each = 4)), x1 = rnorm(120),
                 f = factor(rep(c("a", "b"), 60)))
dc$y <- rnorm(120)
cat("cs(0 + f:x1 | g) cmc TRUE:",
    sh(q(frm(bf(y ~ 1 + cs(0 + f:x1 | g)), data = dc, dry_run = "frame"))),
    "\n")
cat("cs(0 + f:x1 | g) cmc FALSE:",
    sh(q(frm(bf(y ~ 1 + cs(0 + f:x1 | g), cmc = FALSE), data = dc,
             dry_run = "frame"))), "\n")
cat("DONE\n")
