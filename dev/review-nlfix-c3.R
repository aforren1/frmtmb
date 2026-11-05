# Item 3: cost, in one process, minimum of repeats. FRM_ARM=before gives the base for the
# models both builds fit. Seed 6401, the review's. Machine shared; read the minima.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
tm <- function(f, k = 5) { v <- numeric(k); for (i in seq_len(k)) v[i] <- system.time(f())[["elapsed"]]; min(v) }
set.seed(6401)
mk <- function(n, ng) {
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(ng, n, TRUE)))
  d$y <- 2 * exp(0.3 * d$x) + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5)
  d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
  d
}
q <- function(expr) suppressWarnings(suppressMessages(expr))
st <- list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.1))
d1 <- mk(2e5, 200); d2 <- mk(2e4, 200); d3 <- mk(3e4, 30); d4 <- mk(2000, 50)
cases <- list(
  "plain 2e5" = function() q(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d1, start = st)),
  "no start 2e5" = function() q(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), data = d1)),
  "RE 2e4" = function() q(frm(bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE), data = d2, start = st)),
  "plain 2e3" = function() q(frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d4, start = st)))
for (nm in names(cases)) cat(sprintf("%-16s total fit, min of 5: %.3f s\n", nm, tm(cases[[nm]])))
if (arm == "before") quit("no")
cat("\ncumulative 3e4, body against the linear spelling on the same build; no flag\n")
tn <- tm(function() q(frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d3)))
tl <- tm(function() q(frm(yo ~ 0 + x, family = cumulative(), data = d3)))
cat(sprintf("body %.3f s, linear %.3f s, ratio %.3f\n", tn, tl, tn / tl))
fit <- q(frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d3))
nll <- ns$build_objective(fit$frame); tpl <- ns$make_start(fit$frame, NULL)
block <- function(f) { k <- 1L; repeat { t0 <- proc.time()[["elapsed"]]; for (i in seq_len(k)) f(); dt <- proc.time()[["elapsed"]] - t0; if (dt > 0.5) break; k <- k * 2L }; dt / k }
g <- min(replicate(3, block(function() ns$nl_check_body_read(fit$frame, nll, tpl))))
s <- min(replicate(3, block(function() ns$nl_saddle_find(fit))))
dm <- min(replicate(3, block(function() ns$nl_domain_start(fit$frame, tpl, integer(0), FALSE))))
cat(sprintf("guard %.4f domain %.4f saddle look %.4f, all / body fit = %.3f\n", g, dm, s, (g + dm + s) / tn))
cat("\nthe saddle, 3e4 rows, decomposed\n")
fa <- bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
t_all <- system.time(ra <- q(frm(fa, family = cumulative(), data = d3)))[["elapsed"]]
cat(sprintf("frm() default start: %.2f s, ll %.4f, beta %s\n", t_all, as.numeric(logLik(ra)), paste(format(ra$estimates$beta, digits = 4), collapse = ",")))
sad <- list(seen = c("b1", "b2"), lp = fit$frame$linpreds[[1]])
ff <- q(frm(fa, family = cumulative(), data = d3, start = list(beta = c("b1_(Intercept)" = 0, "b2_(Intercept)" = 0)), control = frmtmb_control(check_se = "ignore")))
bd <- Filter(function(q) !is.null(q[["nl_body"]]), ra$frame$linpreds); stt <- ns$nl_saddle_starts(ra$frame, lapply(bd, function(q) list(lp = q, seen = c("b1", "b2"))), NULL)
for (i in seq_along(stt)) {
  a <- system.time(r1 <- q(frm(fa, family = cumulative(), data = d3, start = stt[[i]])))[["elapsed"]]
  b <- system.time(r2 <- q(frm(fa, family = cumulative(), data = d3, start = stt[[i]], se = FALSE)))[["elapsed"]]
  cat(sprintf("restart %d (%s): with se %.2f s, ll %.4f conv %d iter %s | se = FALSE %.2f s, ll %.4f\n", i, paste(format(unlist(stt[[i]]), digits = 3), collapse = ","), a, as.numeric(logLik(r1)), r1$opt$convergence, r1$opt$iterations, b, as.numeric(logLik(r2))))
}
cat("\nan ordinary flagged fit: what the restarts cost, and what they did\n")
set.seed(6402)
dr <- mk(3000, 30); dr$c <- rpois(3000, exp(0.3 + 0.2 * dr$x))
fl <- list(
  "a + b, ridge" = function(s) q(frm(bf(y ~ a + b, a + b ~ 1, nl = TRUE), data = dr, start = s)),
  "RE sd at 0" = function(s) q(frm(bf(y ~ a + b * x, a ~ 1 + (1 | g), b ~ 1, nl = TRUE), data = transform(dr, y = rnorm(3000)), start = s)),
  "shared mixture" = function(s) q(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE), family = mixture(gaussian(), gaussian()), data = dr, start = s)),
  "b0 + b1 x cumul" = function(s) q(frm(bf(yo ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = cumulative(), data = dr, start = s)))
for (nm in names(fl)) {
  t1 <- tm(function() fl[[nm]](NULL), 2)
  # the same fit with start = its default values, which takes no restart
  fit1 <- fl[[nm]](NULL)
  sv <- list(beta = setNames(rep(0, length(fit1$estimates$beta)), names(fit1$estimates$beta)))
  t2 <- tryCatch(tm(function() fl[[nm]](sv), 2), error = function(e) NA)
  cat(sprintf("%-18s default %.2f s, given the zero start %.2f s, ratio %.1f\n", nm, t1, t2, t1 / t2))
}
