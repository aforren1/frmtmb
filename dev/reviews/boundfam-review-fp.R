# Review of lane boundfam: false positives of se_edge_limit() on fits
# with a real, finite maximum far out on the link, and the units probe.
#
#   Rscript dev/reviews/boundfam-review-fp.R <lane|base> [case ...]
#
# Prints per case: conditions, estimates, SEs, se_lost; and a reference
# fit (glm, MASS::glm.nb, glmmTMB) where one exists.
arm <- commandArgs(trailingOnly = TRUE)[1]
cases <- commandArgs(trailingOnly = TRUE)[-1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")

conds <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, paste0("[", class(x)[1], "] ", conditionMessage(x)))
    invokeRestart("muffleMessage")
  })
  list(value = v, w = w, m = m)
}
show <- function(name, r) {
  fit <- r$value
  cat("\n=== ", name, "\n", sep = "")
  for (x in r$w) cat("  W: ", substr(gsub("\n", " ", x), 1, 400), "\n")
  for (x in r$m) cat("  M: ", substr(gsub("\n", " ", x), 1, 400), "\n")
  if (!inherits(fit, "frmtmb_fit")) { cat("  no fit\n"); return(invisible()) }
  est <- fit$opt$par
  V <- suppressWarnings(vcov(fit, full = TRUE))
  se <- suppressWarnings(sqrt(diag(V)))
  nm <- names(se)
  lost <- frmtmb:::sdr_of(fit)$se_lost
  cat(sprintf("  code=%s obj=%.8f\n", fit$opt$convergence, fit$opt$objective))
  for (i in seq_along(se)) {
    cat(sprintf("  %-24s est=%12.6g se=%12.6g %s\n", nm[i], est[i], se[i],
                if (nm[i] %in% names(lost)) lost[[nm[i]]] else ""))
  }
  invisible()
}
ref <- function(name, m) {
  cat("  REF ", name, ": ", sep = "")
  if (inherits(m, "try-error")) { cat("error\n"); return(invisible()) }
  s <- tryCatch(summary(m)$coefficients, error = function(e) NULL)
  if (is.list(s) && !is.matrix(s)) s <- s$cond
  ll <- tryCatch(as.numeric(logLik(m)), error = function(e) NA)
  cat(sprintf("logLik=%.8f\n", ll))
  if (is.matrix(s)) {
    for (i in seq_len(nrow(s))) {
      cat(sprintf("    %-22s %12.6g se=%12.6g\n", rownames(s)[i], s[i, 1], s[i, 2]))
    }
  }
  invisible()
}
ok <- function(c) !length(cases) || c %in% cases
F <- function(...) conds(frm(...))

# 1. rare-event bernoulli, finite MLE far out on the logit
if (ok("bern_rare")) for (s in 1:3) {
  set.seed(s)
  n <- 1e5
  d <- data.frame(x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-10 + 0.5 * d$x))
  cat("\nbern_rare seed", s, "events", sum(d$y), "\n")
  show(paste0("bern_rare s", s), F(y ~ x, family = bernoulli(), data = d))
  ref("glm", glm(y ~ x, binomial, d))
}
# intercept only, 1 and 2 events in 2e5 (intercept -12.2, -11.5)
if (ok("bern_one")) for (k in 1:2) {
  d <- data.frame(y = c(rep(1L, k), rep(0L, 2e5 - k)))
  show(paste0("bern_one k", k), F(y ~ 1, family = bernoulli(), data = d))
  ref("glm", glm(y ~ 1, binomial, d))
}
# 2. rare-event poisson
if (ok("pois_rare")) for (s in 1:3) {
  set.seed(s)
  n <- 2e4
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(-9 + 0.3 * d$x))
  cat("\npois_rare seed", s, "events", sum(d$y), "\n")
  if (sum(d$y) == 0) d$y[1] <- 1L
  show(paste0("pois_rare s", s), F(y ~ x, family = poisson(), data = d))
  ref("glm", glm(y ~ x, poisson, d))
}
# 3. quasi-separation that is not complete: one overlap
if (ok("quasi")) for (s in 1:3) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- as.integer(d$x > 0)
  i <- which.min(abs(d$x - 0.3)); d$y[i] <- 0L
  j <- which.min(abs(d$x + 0.3)); d$y[j] <- 1L
  show(paste0("quasi s", s), F(y ~ x, family = bernoulli(), data = d))
  ref("glm", glm(y ~ x, binomial, d))
}
# 4. a factor level with 3 rows and 1 event (finite), and 0 events
if (ok("fac3")) for (ev in c(1L, 0L)) {
  set.seed(4)
  n <- 3000
  f <- factor(c(rep("a", n - 3), rep("b", 3)))
  d <- data.frame(f = f, x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-4))
  d$y[f == "b"] <- c(rep(1L, ev), rep(0L, 3 - ev))
  show(paste0("fac3 bern ev", ev), F(y ~ x + f, family = bernoulli(), data = d))
  ref("glm", glm(y ~ x + f, binomial, d))
  d$y <- rpois(n, exp(-2))
  d$y[f == "b"] <- c(rep(1L, ev), rep(0L, 3 - ev))
  show(paste0("fac3 pois ev", ev), F(y ~ x + f, family = poisson(), data = d))
  ref("glm", glm(y ~ x + f, poisson, d))
}
# 5. negbinomial with a large but finite shape, large n
if (ok("nb_large")) for (sh in c(50, 300)) for (s in 1:2) {
  set.seed(s)
  n <- 2e4
  d <- data.frame(x = rnorm(n))
  d$y <- MASS::rnegbin(n, exp(1 + 0.3 * d$x), sh)
  show(paste0("nb shape", sh, " s", s), F(y ~ x, family = negbinomial(), data = d))
  ref("glm.nb", try(MASS::glm.nb(y ~ x, d), silent = TRUE))
  g <- try(MASS::glm.nb(y ~ x, d), silent = TRUE)
  if (!inherits(g, "try-error")) cat(sprintf("    theta=%.6g se=%.6g log=%.6g\n", g$theta, g$SE.theta, log(g$theta)))
  ref("glm pois", glm(y ~ x, poisson, d))
}
# 6. beta with phi large but finite
if (ok("beta_large")) for (ph in c(500, 5000)) {
  set.seed(6)
  n <- 2000
  d <- data.frame(x = rnorm(n))
  m <- plogis(0.2 + 0.3 * d$x)
  d$y <- rbeta(n, m * ph, (1 - m) * ph)
  show(paste0("beta phi", ph), F(y ~ x, family = Beta(), data = d))
  ref("glmmTMB", try(glmmTMB::glmmTMB(y ~ x, d, family = glmmTMB::beta_family()), silent = TRUE))
}
# 7. student nu around 50 to 200, large n
if (ok("student_nu")) for (nu in c(50, 200)) for (s in 1:2) {
  set.seed(s)
  n <- 2e4
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rt(n, nu)
  show(paste0("student nu", nu, " s", s), F(y ~ x, family = student(), data = d))
  ref("glmmTMB t", try(glmmTMB::glmmTMB(y ~ x, d, family = glmmTMB::t_family()), silent = TRUE))
  g <- try(glmmTMB::glmmTMB(y ~ x, d, family = glmmTMB::t_family()), silent = TRUE)
  if (!inherits(g, "try-error")) cat(sprintf("    glmmTMB log nu=%.6g\n", log(glmmTMB::family_params(g))))
}
# 8. zero-inflation logit near -8 with a real maximum
if (ok("zi8")) for (s in 1:3) {
  set.seed(s)
  n <- 5e4
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, plogis(-8))
  d$y <- ifelse(z == 1, 0L, rpois(n, exp(2.2 + 0.1 * d$x)))
  cat("\nzi8 seed", s, "structural zeros", sum(z), "zeros", sum(d$y == 0), "\n")
  show(paste0("zi8 s", s), F(y ~ x, family = zero_inflated_poisson(), data = d))
  ref("glmmTMB zip", try(glmmTMB::glmmTMB(y ~ x, d, family = poisson, ziformula = ~1), silent = TRUE))
  g <- try(glmmTMB::glmmTMB(y ~ x, d, family = poisson, ziformula = ~1), silent = TRUE)
  if (!inherits(g, "try-error")) print(summary(g)$coefficients$zi)
}
# 9. units: covariate at extreme scales, autoscale on (default) and off,
# with and without a smooth (which sets par_units without autoscale)
if (ok("units")) for (sc in c(1e-6, 1e6)) for (as in list(NULL, FALSE)) {
  set.seed(9)
  n <- 300
  d <- data.frame(x = rnorm(n), z = runif(n))
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + sin(2 * pi * d$z)))
  d$xs <- d$x * sc
  ctl <- frmtmb_control(autoscale = as)
  lab <- paste0("units sc=", sc, " autoscale=", format(as))
  show(paste(lab, "y ~ xs"), F(y ~ xs, family = poisson(), data = d, control = ctl))
  show(paste(lab, "y ~ xs + s(z)"), F(y ~ xs + s(z), family = poisson(), data = d, control = ctl))
  ref("glm", glm(y ~ xs, poisson, d))
}
# 10. barely interior negbinomial on poisson data: how often fires,
# and the poisson logLik next to the negbinomial one where it fires
if (ok("nb_pois")) for (s in 1:30) {
  set.seed(s)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- rpois(n, exp(1 + 0.4 * d$x))
  r <- F(y ~ x, family = negbinomial(), data = d)
  fit <- r$value
  lost <- frmtmb:::sdr_of(fit)$se_lost
  gp <- glm(y ~ x, poisson, d)
  g <- suppressWarnings(try(MASS::glm.nb(y ~ x, d), silent = TRUE))
  cat(sprintf("nb_pois s%02d logshape=%8.4f lost=%-12s nW=%d nM=%d ll_frm=%.8f ll_pois=%.8f ll_nb=%s theta_nb=%s\n",
              s, fit$opt$par[length(fit$opt$par)], paste(lost, collapse = ","), length(r$w), length(r$m),
              -fit$opt$objective, as.numeric(logLik(gp)),
              if (inherits(g, "try-error")) "err" else sprintf("%.8f", as.numeric(logLik(g))),
              if (inherits(g, "try-error")) "err" else sprintf("%.4g", g$theta)))
}
# 11. student nu: n = 200 normal data, how often fires, gaussian ll
if (ok("st_norm")) for (s in 1:20) {
  set.seed(s)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rt(n, 60)
  r <- F(y ~ x, family = student(), data = d)
  fit <- r$value
  lost <- frmtmb:::sdr_of(fit)$se_lost
  gg <- glm(y ~ x, gaussian, d)
  se <- suppressWarnings(sqrt(diag(vcov(fit, full = TRUE))))
  k <- grep("^nu", names(se))[1]
  cat(sprintf("st_norm s%02d %s=%8.4f se=%9.4g lost=%-12s nW=%d nM=%d ll_frm=%.8f ll_gauss=%.8f\n",
              s, names(se)[k], fit$opt$par[k], se[k], paste(lost, collapse = ","), length(r$w), length(r$m),
              -fit$opt$objective, as.numeric(logLik(gg))))
}
