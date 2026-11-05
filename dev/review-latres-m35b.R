# Short check, item 5: M3 and M5. Seed 5 at sderr = 0; the message's marginal sd;
# diagnose()$singular; the limit under ar -> -1 and a random walk; whether the SE
# check standing aside hides another problem; nlf(zi ~ ...) by a Laplace identity.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
mkd <- function(seed, dgp = "ri", fam = "pois", G = 10, K = 6) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- rnorm(n)
  re <- switch(dgp,
    ri = rep(rnorm(G, 0, .5), each = K),
    none = 0,
    walk = as.vector(replicate(G, cumsum(rnorm(K, 0, .4)))),
    alt = as.vector(replicate(G, rnorm(1, 0, .6) * rep(c(1, -1), length.out = K))))
  eta <- .5 + .3 * d$x + re
  d$y <- rpois(n, exp(eta))
  d$bin <- rbinom(n, 1, plogis(eta))
  d
}
capt <- function(expr) {
  msgs <- character(0)
  wr <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }, warning = function(w) {
    wr <<- c(wr, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(v = v, m = msgs, w = wr)
}
show_fit <- function(lab, r) {
  f <- r$v
  cat("\n---", lab, "\n")
  if (inherits(f, "error")) { cat("ERR", conditionMessage(f), "\n"); return(invisible()) }
  th <- f$estimates$theta
  cat("code", f$opt$convergence, " ll", as.numeric(logLik(f)), " theta", th, "\n")
  for (m in r$m) cat("MSG:", substr(m, 1, 400), "\n")
  for (m in r$w) cat("WARN:", substr(m, 1, 300), "\n")
  print(summary(f)$cor_pars)
  print(diagnose(f)$singular)
  print(confint_varcorr(f))
  cat("fixed SE", sqrt(diag(as.matrix(vcov(f)))), "\n")
}
show_fit("RI, seed 5, sderr = 0", capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = mkd(5))))
fx <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = mkd(1)))
show_fit("RI, seed 1, ar = 1", fx)
f1 <- fx$v
th <- f1$estimates$theta
cat("marginal sd by hand:", exp(th[1]) * sqrt(1 + th[2]^2), " (1 | g) sd:", {
  fr <- frm(y ~ x + (1 | g), family = poisson(), data = mkd(1))
  sqrt(as.numeric(VarCorr(fr)$g$sd[1])^2)
}, " cosy sderr:", exp(frm(y ~ x + cosy(t, g), family = poisson(), data = mkd(1))$estimates$theta[1]), "\n")

cat("\n===== is the sentence true elsewhere? search seeds for a limit under other shapes\n")
tally <- function(dgp, fam = "pois", seeds = 1:40, G = 10, K = 6, term = "ar(t, g, cov = TRUE)") {
  hits <- 0L
  first <- NULL
  for (s in seeds) {
    d <- mkd(s, dgp, G = G, K = K)
    resp <- if (fam == "pois") "y" else "bin"
    r <- capt(frm(as.formula(paste(resp, "~ x +", term)), family = if (fam == "pois") poisson() else bernoulli(), data = d))
    f <- r$v
    if (inherits(f, "error")) next
    bkx <- fi("frame_latres")(f$frame)[[1]]; a <- bkx$aux
    lim <- fi("latres_block_limit")(f, bkx)$at
    if (any(lim)) {
      hits <- hits + 1L
      tt <- f$estimates$theta
      ar <- tt[2] / sqrt(1 + tt[2]^2)
      read <- any(grepl("read as cosy", r$m))
      ll_c <- capt(frm(as.formula(paste(resp, "~ x + (1 | g)")), family = if (fam == "pois") poisson() else bernoulli(), data = d))$v
      gapc <- as.numeric(logLik(f)) - as.numeric(logLik(ll_c))
      cat(sprintf("  %-5s seed %2d ar %7.4f sderr %.4g  'read as cosy' %-5s  ll - ll(1|g) %.3f  code %d\n", dgp, s, ar, exp(tt[1]), read, gapc, f$opt$convergence))
    }
  }
  cat(dgp, fam, ": at a limit", hits, "of", length(seeds), "\n")
}
if (FALSE) tally("walk", K = 12)
tally("alt")
if (FALSE) tally("none")
if (FALSE) tally("walk", "bern", K = 20, G = 15, seeds = 1:15)

cat("\n===== the SE check standing aside: a fit at a limit with another lost parameter\n")
# a design alias beside the limit: x2 = x, and a beta with no curvature
d <- mkd(1)
d$x2 <- d$x * 1
d$lone <- as.numeric(d$g == "1" & d$t == 1)
d$sep <- as.numeric(d$y == 0)
r <- capt(frm(y ~ x + sep + ar(t, g, cov = TRUE), family = poisson(), data = d))
show_fit("limit plus a separated slope", r)
r <- capt(frm(y ~ x + x2 + ar(t, g, cov = TRUE), family = poisson(), data = d))
show_fit("limit plus a design alias", r)
d2 <- mkd(1)
d2$nb <- d2$y
r <- capt(frm(nb ~ x + ar(t, g, cov = TRUE), family = negbinomial(), data = d2))
show_fit("limit with shape to infinity (no overdispersion in the rows)", r)
r <- capt(frm(y ~ x + (x | g) + ar(t, g, cov = TRUE), family = poisson(), data = d2))
show_fit("limit beside a random slope at its own limit", r)
cat("\n===== with check_se\n")
for (lev in c("ignore", "warning", "message", "error")) {
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = mkd(1), control = frmtmb_control(check_se = lev)))
  cat(lev, ": err", if (inherits(r$v, "error")) substr(conditionMessage(r$v), 1, 100) else "none", " msgs", length(r$m), " warns", length(r$w), "\n")
}
