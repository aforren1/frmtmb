# Reviewer, item 1: the merge resolutions on ordinal mixtures.
#  (a) default_prior() Intercept/delta rows against brms 2.23.0 on
#      every mixture shape that both take;
#  (b) per-threshold coef priors on structures that have no such row;
#  (c) the log prior density of per-threshold priors on mixture
#      components against a hand computation (an independent shape);
#  (d) confint() / vcov(full = TRUE) names: equal, unique, and each
#      resolves through confint(parm = ) to its own row.
# Data seed 20261006.
#   Rscript dev/relrev-mix.R > dev/relrev-log/mix.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261006)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
d$yh <- ifelse(runif(n) < 0.15, 0L, d$y)
quiet <- function(expr) {
  tryCatch(suppressMessages(suppressWarnings(expr)),
           error = function(e) structure(conditionMessage(e), class = "err"))
}
rows <- function(p) {
  p <- as.data.frame(p)
  keep <- p$class %in% c("Intercept", "delta")
  sort(paste(p$class, p$coef, p$group, p$dpar, sep = "|")[keep])
}
shapes <- list(
  none_cum_sratio = list(f = quote(mixture(cumulative(), sratio())),
                         b = quote(brms::mixture(brms::cumulative(), brms::sratio()))),
  none_cum_cum = list(f = quote(mixture(cumulative(), cumulative())),
                      b = quote(brms::mixture(brms::cumulative(), brms::cumulative()))),
  none_sratio_acat = list(f = quote(mixture(sratio(), acat())),
                          b = quote(brms::mixture(brms::sratio(), brms::acat()))),
  none_three = list(f = quote(mixture(cumulative(), sratio(), cratio())),
                    b = quote(brms::mixture(brms::cumulative(), brms::sratio(), brms::cratio()))),
  mu_cum_sratio = list(f = quote(mixture(cumulative(), sratio(), order = "mu")),
                       b = quote(brms::mixture(brms::cumulative(), brms::sratio(), order = "mu"))),
  mu_sratio_acat = list(f = quote(mixture(sratio(), acat(), order = "mu")),
                        b = quote(brms::mixture(brms::sratio(), brms::acat(), order = "mu"))),
  none_equi_cum = list(f = quote(mixture(cumulative(threshold = "equidistant"), cumulative())),
                       b = quote(brms::mixture(brms::cumulative(threshold = "equidistant"), brms::cumulative()))),
  none_stz_sratio = list(f = quote(mixture(cumulative(threshold = "sum_to_zero"), sratio())),
                         b = quote(brms::mixture(brms::cumulative(threshold = "sum_to_zero"), brms::sratio()))),
  mu_stz_cum = list(f = quote(mixture(cumulative(threshold = "sum_to_zero"), cumulative(), order = "mu")),
                    b = quote(brms::mixture(brms::cumulative(threshold = "sum_to_zero"), brms::cumulative(), order = "mu"))),
  none_hurdle_cum = list(f = quote(mixture(hurdle_cumulative(), hurdle_cumulative())),
                         b = quote(brms::mixture(brms::hurdle_cumulative(), brms::hurdle_cumulative())), y = "yh")
)
cat("\n== (a) default_prior() threshold rows, frmtmb against brms\n")
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  yv <- s$y %||% "y"
  fo <- stats::as.formula(paste(yv, "~ x"))
  fr <- quiet(rows(default_prior(bf(fo), data = d, family = eval(s$f))))
  br <- quiet(rows(brms::default_prior(brms::bf(fo), data = d, family = eval(s$b))))
  same <- !inherits(fr, "err") && !inherits(br, "err") && identical(fr, br)
  cat(sprintf("%-18s %s  frmtmb %s | brms %s\n", nm, if (same) "SAME" else "DIFF",
              if (inherits(fr, "err")) paste("ERR", substr(fr, 1, 80)) else paste(length(fr), "rows"),
              if (inherits(br, "err")) paste("ERR", substr(br, 1, 80)) else paste(length(br), "rows")))
  if (!same && !inherits(fr, "err") && !inherits(br, "err")) {
    cat("   only frmtmb:", setdiff(fr, br), "\n   only brms:  ", setdiff(br, fr), "\n")
  }
}
cat("\n== (b) coef priors where brms has no such row\n")
pb <- function(tag, fam, pr, y = "y") {
  r <- quiet(frm(bf(stats::as.formula(paste(y, "~ x"))), family = fam, data = d, prior = pr))
  cat(sprintf("%-46s %s\n", tag, if (inherits(r, "err")) paste("REFUSED:", substr(r, 1, 110)) else "ACCEPTED"))
}
pb("none equi comp 1, coef 2 dpar mu1", mixture(cumulative(threshold = "equidistant"), cumulative()),
   set_prior("normal(0, 1)", class = "Intercept", coef = "2", dpar = "mu1"))
pb("none stz comp 1, coef 2 dpar mu1", mixture(cumulative(threshold = "sum_to_zero"), sratio()),
   set_prior("normal(0, 1)", class = "Intercept", coef = "2", dpar = "mu1"))
pb("none, coef 4 of 3 dpar mu2", mixture(cumulative(), sratio()),
   set_prior("normal(0, 1)", class = "Intercept", coef = "4", dpar = "mu2"))
pb("none, coef 2 with no dpar", mixture(cumulative(), sratio()),
   set_prior("normal(0, 1)", class = "Intercept", coef = "2"))
pb("mu, coef 2 dpar mu1", mixture(cumulative(), sratio(), order = "mu"),
   set_prior("normal(0, 1)", class = "Intercept", coef = "2", dpar = "mu1"))
pb("none, coef 2 dpar mu3 (no such comp)", mixture(cumulative(), sratio()),
   set_prior("normal(0, 1)", class = "Intercept", coef = "2", dpar = "mu3"))
cat("brms on the same rows (stancode with the prior):\n")
bb <- function(tag, fam, pr) {
  r <- quiet(brms::stancode(brms::bf(y ~ x), data = d, family = fam, prior = pr))
  cat(sprintf("  brms %-41s %s\n", tag, if (inherits(r, "err")) paste("REFUSED:", substr(gsub("\n", " ", r), 1, 100)) else "ACCEPTED"))
}
bb("none, coef 2 with no dpar", brms::mixture(brms::cumulative(), brms::sratio()),
   brms::set_prior("normal(0, 1)", class = "Intercept", coef = "2"))
bb("mu, coef 2 dpar mu1", brms::mixture(brms::cumulative(), brms::sratio(), order = "mu"),
   brms::set_prior("normal(0, 1)", class = "Intercept", coef = "2", dpar = "mu1"))
bb("none, coef 4 of 3 dpar mu2", brms::mixture(brms::cumulative(), brms::sratio()),
   brms::set_prior("normal(0, 1)", class = "Intercept", coef = "4", dpar = "mu2"))

cat("\n== (c) log prior density against a hand computation\n")
lp_of <- function(fit) {
  ent <- frmtmb:::resolve_prior_input(list(frame = fit$frame, spec = fit$spec),
                                      fit$prior)$entries
  -frmtmb:::neg_log_prior_fn(ent)(fit$estimates)
}
# order none, sratio first and cumulative second: a coef prior on each,
# plus a class prior on component 1 that the coef row overrides at its
# threshold. brms centers each component's design: its Intercept_mu<k>
# is b_mu<k>_Intercept - mean(x) * b_mu<k>_x
pr <- set_prior("normal(0.2, 0.7)", class = "Intercept", coef = "1", dpar = "mu2") +
  set_prior("normal(1, 0.5)", class = "Intercept", coef = "3", dpar = "mu1") +
  set_prior("student_t(3, 0, 2.5)", class = "Intercept", dpar = "mu1")
fit <- quiet(frm(bf(y ~ x), family = mixture(sratio(), cumulative()), data = d, prior = pr))
if (inherits(fit, "err")) cat("fit ERR", fit, "\n") else {
  e <- fit$estimates
  cat("estimate names:", names(e), "\n")
  t1 <- e$tau_raw1                                   # sratio: thresholds as held
  r2 <- e$tau_raw2; t2 <- c(r2[1], r2[1] + cumsum(exp(r2[-1])))
  b1 <- e$beta[["mu1_x"]]; b2 <- e$beta[["mu2_x"]]
  m <- mean(d$x)
  c1 <- t1 - m * b1; c2 <- t2 - m * b2
  dst <- function(x, nu, mu, s) dt((x - mu) / s, nu, log = TRUE) - log(s)
  ref <- sum(dst(c1[1:2], 3, 0, 2.5)) + dnorm(c1[3], 1, 0.5, log = TRUE) +
    dnorm(c2[1], 0.2, 0.7, log = TRUE) + sum(r2[-1])
  got <- lp_of(fit)
  cat(sprintf("none sratio+cum: frmtmb %.15g  hand %.15g  rel %.3g\n", got, ref, abs(got - ref) / abs(ref)))
  # the same hand density without the centering, and without the
  # Jacobian, to show the check can fail
  ref_nc <- sum(dst(t1[1:2], 3, 0, 2.5)) + dnorm(t1[3], 1, 0.5, log = TRUE) +
    dnorm(t2[1], 0.2, 0.7, log = TRUE) + sum(r2[-1])
  cat(sprintf("   uncentered hand %.15g (diff %.3g); no Jacobian diff %.3g\n", ref_nc, got - ref_nc, sum(r2[-1])))
}
# order mu, unordered components only (sratio, acat): no Jacobian, no centering
fit2 <- quiet(frm(bf(y ~ x), family = mixture(sratio(), acat(), order = "mu"), data = d,
                  prior = set_prior("normal(-0.5, 0.4)", class = "Intercept", coef = "2")))
if (inherits(fit2, "err")) cat("fit2 ERR", fit2, "\n") else {
  t <- fit2$estimates$tau_raw
  ref <- dnorm(t[2], -0.5, 0.4, log = TRUE)
  got <- lp_of(fit2)
  cat(sprintf("mu sratio+acat: frmtmb %.15g  hand %.15g  rel %.3g\n", got, ref, abs(got - ref) / abs(ref)))
}
cat("\n== (d) confint() and vcov(full = TRUE) names on mixture shapes\n")
fshapes <- list(
  none_cum_sratio = quote(mixture(cumulative(), sratio())),
  none_cum_cum = quote(mixture(cumulative(), cumulative())),
  none_sratio_acat = quote(mixture(sratio(), acat())),
  none_three = quote(mixture(cumulative(), sratio(), cratio())),
  mu_cum_sratio = quote(mixture(cumulative(), sratio(), order = "mu")),
  mu_sratio_acat = quote(mixture(sratio(), acat(), order = "mu")),
  mu_cum_cum = quote(mixture(cumulative(), cumulative(), order = "mu")),
  none_equi_cum = quote(mixture(cumulative(threshold = "equidistant"), cumulative())),
  none_equi_acat = quote(mixture(acat(threshold = "equidistant"), sratio())),
  none_stz_sratio = quote(mixture(cumulative(threshold = "sum_to_zero"), sratio())),
  mu_stz_cum = quote(mixture(cumulative(threshold = "sum_to_zero"), cumulative(), order = "mu")),
  none_disc = quote(mixture(cumulative(), cumulative())),
  none_hurdle = quote(mixture(hurdle_cumulative(), cumulative())),
  none_cs = quote(mixture(cumulative(), sratio()))
)
for (nm in names(fshapes)) {
  fo <- switch(nm, none_disc = bf(y ~ x, disc1 ~ 0 + z),
               none_hurdle = bf(yh ~ x), none_cs = bf(y ~ cs(z) + x), bf(y ~ x))
  fit <- quiet(frm(fo, family = eval(fshapes[[nm]]), data = d))
  if (inherits(fit, "err")) {
    cat(sprintf("%-16s FIT ERR %s\n", nm, substr(fit, 1, 120)))
    next
  }
  ci <- quiet(confint(fit))
  vc <- quiet(vcov(fit, full = TRUE))
  if (inherits(ci, "err") || inherits(vc, "err")) {
    cat(sprintf("%-16s confint/vcov ERR %s %s\n", nm, ci, vc)); next
  }
  rn <- rownames(ci)
  thr <- rn[grepl("Intercept|delta|tau", rn)]
  eqv <- identical(rownames(vc), rn)
  uniq <- !anyDuplicated(rn)
  # each label resolves through parm = to its own row (est and bounds)
  res <- vapply(rn, function(p) {
    r <- quiet(suppressMessages(confint(fit, parm = p)))
    if (inherits(r, "err")) return("ERR")
    if (nrow(r) != 1L) return(paste0("n=", nrow(r)))
    if (identical(unname(r[1, ]), unname(ci[p, ]))) "ok" else "DIFF"
  }, "")
  cat(sprintf("%-16s rows %2d unique %s vcov-same %s parm %s; thresholds: %s\n",
              nm, length(rn), uniq, eqv,
              if (all(res == "ok")) "all-ok" else paste(names(res)[res != "ok"], res[res != "ok"], collapse = ";"),
              paste(thr, collapse = ", ")))
}
