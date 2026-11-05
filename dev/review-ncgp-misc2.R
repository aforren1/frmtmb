# Items 3, 5: does anything reach gp_fd_step() and move when it is a no-op; readers spelled
# right; kriging through nl; frmtmb.spatial; by=, mvbf, sigma ~ gp() against base numbers.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
orig_step <- ns$gp_fd_step
calls <- 0L
counted <- function(g, fit) { calls <<- calls + 1L; orig_step(g, fit) }
noop <- function(g, fit) g
mk <- function(fam, n, seed, sg = NA) {
  set.seed(seed); x <- sort(runif(n, 0, 8))
  eta <- 1.5 + 1.2 * sin(x)
  d <- data.frame(x = x, z = rnorm(n))
  switch(fam,
    pois = { d$y <- rpois(n, exp(eta)); list(d = d, fml = y ~ gp(x), fam = poisson()) },
    pois6 = { d$y <- rpois(n, exp(eta + 6)); list(d = d, fml = y ~ gp(x), fam = poisson()) },
    gamma = { d$y <- rgamma(n, 20, 20 / exp(eta)); list(d = d, fml = y ~ gp(x), fam = Gamma(link = "log")) },
    lnorm = { d$y <- exp(eta + rnorm(n, sd = 0.01)); list(d = d, fml = y ~ gp(x), fam = lognormal()) },
    ord = { d$y <- ordered(cut(eta + rlogis(n), c(-Inf, 1, 2, 3, Inf), labels = FALSE)); list(d = d, fml = y ~ gp(x), fam = cumulative()) })
}
for (cs in list(list("pois", 60), list("pois", 200), list("pois6", 200), list("gamma", 200), list("lnorm", 200), list("ord", 150))) {
  m <- mk(cs[[1]], cs[[2]], 900 + cs[[2]])
  f <- suppressWarnings(frm(m$fml, m$d, family = m$fam))
  nd <- data.frame(x = c(0.5, 2.1, 4.33, 7.4), z = 0)
  rd <- function() {
    fcopy <- f; fcopy$cache <- new.env(parent = emptyenv())
    list(a = tryCatch(suppressWarnings(fitted(fcopy)), error = function(e) NA), b = tryCatch(suppressWarnings(fitted(fcopy, newdata = nd)), error = function(e) NA))
  }
  calls <- 0L; patch("gp_fd_step", counted); r0 <- rd(); nc <- calls
  patch("gp_fd_step", noop); r1 <- rd(); patch("gp_fd_step", orig_step)
  df <- function(a, b) { a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b)); ok <- is.finite(a) & is.finite(b); if (!any(ok)) return(NA); max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-8)) }
  cat(sprintf("RMV %s n %d conv %d gp_fd_step called %d | no-op changes fitted %.1e newdata %.1e\n", cs[[1]], cs[[2]], f$opt$convergence, nc, df(r1$a, r0$a), df(r1$b, r0$b)))
  flush.console()
}
rr <- function(lab, e) cat(sprintf("RRD %-34s %s\n", lab, tryCatch({ v <- suppressWarnings(suppressMessages(e)); if (is.data.frame(v)) v <- as.matrix(Filter(is.numeric, v)); paste(class(v)[1], paste(utils::head(signif(suppressWarnings(as.numeric(unlist(v))), 5), 6), collapse = " ")) }, error = function(err) paste("ERR", substr(conditionMessage(err), 1, 150)))))
set.seed(8150); n <- 90; dd <- data.frame(x = runif(n, 0, 6)); dd$y <- sin(dd$x) + rnorm(n, sd = 0.3)
fe <- frm(y ~ gp(x), dd)
xb <- mean(dd$x)
em <- as.data.frame(emmeans::emtrends(fe, ~ 1, var = "x"))
grid <- data.frame(x = c(xb, xb + 0.5))
de <- as.data.frame(frm_curve_deriv(fe, var = "x", newdata = grid))
cat(sprintf("REM emtrends %.5f se %.5f | frm_curve_deriv %.5f se %.5f at x %.4f\n", em$x.trend, em$SE, de$.estimate[1], de$.se[1], xb))
mf <- as.data.frame(marginaleffects::slopes(fe, newdata = data.frame(x = xb), variables = "x"))
cat(sprintf("REM marginaleffects slopes %.5f se %.5f\n", mf$estimate[1], mf$std.error[1]))
rr("frm_curve x", frm_curve(fe, newdata = data.frame(x = c(1, 2, 3))))
rr("frm_curve_feature peak", frm_curve_feature(fe, var = "x", feature = "max", newdata = data.frame(x = seq(0.5, 3, 0.25))))
# through an nl body, fit at low noise
dn <- data.frame(x = runif(50, 0, 6), z = rnorm(50)); dn$y <- 1 + sin(dn$x) * exp(0.2 * dn$z) + rnorm(50, sd = .05)
fn <- suppressWarnings(frm(bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE), dn))
rr("nl fitted", fitted(fn)); rr("nl fitted new", fitted(fn, newdata = dn[1:3, ])); rr("nl ranef", ranef(fn))
rr("nl frm_linpred nlpar", frm_linpred(fn, nlpar = "a", newdata = dn[1:3, ], se.fit = TRUE))
if (requireNamespace("frmtmb.spatial", quietly = TRUE)) {
  library(frmtmb.spatial)
  cat("spatial lib:", dirname(system.file(package = "frmtmb.spatial")), "\n")
  set.seed(77)
  dsp <- data.frame(s1 = runif(40), s2 = runif(40), reg = factor(rep(1:10, 4))); dsp$y <- sin(3 * dsp$s1) + rnorm(40, sd = 0.2)
  W <- matrix(0, 10, 10); for (i in 1:9) W[i, i + 1] <- W[i + 1, i] <- 1; dimnames(W) <- list(1:10, 1:10)
  fs <- tryCatch(suppressWarnings(frm(y ~ gp(s1, s2) + car(W, gr = reg, type = "esicar"), dsp, data2 = list(W = W))), error = function(e) e)
  if (inherits(fs, "error")) cat("spatial err", conditionMessage(fs), "\n") else {
    cat(sprintf("RSP ll %.4f conv %d order %s\n", logLik(fs), fs$opt$convergence, paste(names(fs$obj$env$parameters), collapse = ",")))
    rr("spatial fitted", fitted(fs)); rr("spatial ranef", ranef(fs))
    bk <- fs$frame$re_blocks
    # the field through expand_b against a block by hand: an esicar's b mean removed
    rr("spatial coef_b", ns$coef_b(fs))
  }
}
cat("DONE\n")
