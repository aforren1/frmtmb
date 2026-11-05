# Area 2: extension families and the guard's cost on them. Seed 6301.
source("dev/review-nlfix-common.R")
suppressPackageStartupMessages({ library(frmtmb.eam); library(frmtmb.learn); library(frmtmb.survival); library(frmtmb.latent) })
ns <- asNamespace("frmtmb")
cat("eam from", dirname(find.package("frmtmb.eam")), "learn from", dirname(find.package("frmtmb.learn")), "\n")
ll <- function(r) if (iserr(r)) paste("ERR", substr(r$msg, 1, 230)) else paste0(format(as.numeric(logLik(r)), digits = 10), " conv=", r$opt$convergence, " w=", length(attr(r, "w")))
guard_time <- function(form, fam, data, ...) {
  # the guard's and the domain check's share of this fit, in one process
  dr <- tryCatch(suppressWarnings(suppressMessages(frm(form, family = fam, data = data, dry_run = "objective", ...))), error = function(e) e)
  t_fit <- system.time(fit <- suppressWarnings(suppressMessages(frm(form, family = fam, data = data, ...))))[["elapsed"]]
  fr <- fit$frame
  nll <- ns$build_objective(fr)
  tpl <- ns$make_start(fr, list(...)$start)
  t1 <- system.time(ns$nl_check_body_read(fr, nll, tpl))[["elapsed"]]
  t2 <- system.time(ns$nl_domain_start(fr, tpl, integer(0), FALSE))[["elapsed"]]
  cat(sprintf("   guard %.3f s, domain %.3f s, fit %.3f s, share %.3f\n", t1, t2, t_fit, (t1 + t2) / t_fit))
}
cat("\n--- wiener_gng\n")
set.seed(6301)
d <- wiener_gng_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.25, deadline = 1.5)
d$x <- rnorm(nrow(d))
fam <- wiener_gng(deadline = 1.5)
lin <- tt(frm(bf(rt | dec(responded) ~ x, bias = 0.5), family = fam, data = d))
nl <- tt(frm(bf(rt | dec(responded) ~ b0 + b1 * x, b0 + b1 ~ 1, bias = 0.5, nl = TRUE), family = fam, data = d))
int <- tt(frm(bf(rt | dec(responded) ~ 1, bias = 0.5), family = fam, data = d))
cat("linear", ll(lin), "\nnl", ll(nl), "\nintercept", ll(int), "\n")
if (!iserr(nl)) print(tt(frm_compat("nl", "wiener_gng")))
if (!iserr(nl)) guard_time(bf(rt | dec(responded) ~ b0 + b1 * x, b0 + b1 ~ 1, bias = 0.5, nl = TRUE), fam, d)
cat("\n--- rdm(3), lba(3): several locations, no brms form\n")
dr <- rdm_simulate(150, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
dr$x <- rnorm(nrow(dr))
r <- tt(frm(bf(rt | vint(choice) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = rdm(3), data = dr)); cat("rdm nl =", ll(r), "\n")
r <- tt(frm(bf(rt | vint(choice) ~ 1) + nlf(v2 ~ b0 + b1 * x) + lf(b0 ~ 1) + lf(b1 ~ 1), family = rdm(3), data = dr)); cat("rdm nlf(v2) =", ll(r), "\n")
l <- tt(frm(bf(rt | vint(choice) ~ 1, v2 ~ x), family = rdm(3), data = dr)); cat("rdm linear v2 =", ll(l), "\n")
dl <- lba_simulate(150, v = c(3, 2, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
dl$x <- rnorm(nrow(dl))
r <- tt(frm(bf(rt | vint(choice) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = lba(3), data = dl)); cat("lba nl =", ll(r), "\n")
r <- tt(frm(bf(rt | vint(choice) ~ 1) + nlf(v2 ~ b0 + b1 * x) + lf(b0 ~ 1) + lf(b1 ~ 1), family = lba(3), data = dl)); cat("lba nlf(v2) =", ll(r), "\n")
l <- tt(frm(bf(rt | vint(choice) ~ 1, v2 ~ x), family = lba(3), data = dl)); cat("lba linear v2 =", ll(l), "\n")
cat("\n--- gddm\n")
set.seed(6302)
dg <- gddm_simulate(120L, mu = 2, bs = 2.5, ndt = 0.25, control = gddm_control(t_max = 2))
dg$cond <- rep(1:2, length.out = nrow(dg))
print(names(dg))
cat("\n--- learn: bandit2arm_delta\n")
dt <- frm_task_design("bandit2arm", n_subject = 6L, n_trial = 60L, seed = 11L)
dt$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial), dt, pars = list(alpha = 0.4, tau = 3), seed = 11L)[[1L]]$choice
sub <- as.integer(factor(dt$id)); zs <- c(-1, -0.5, -0.2, 0.2, 0.5, 1); dt$z <- zs[sub]
fm <- bandit2arm_delta(subject = id, trial = trial)
int <- tt(frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1), family = fm, data = dt)); cat("intercept", ll(int), "\n")
lin <- tt(frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ z), family = fm, data = dt)); cat("linear tau ~ z", ll(lin), "\n")
nl1 <- tt(frm(bf(choice | reward(pay1, pay2) ~ 1) + nlf(tau ~ t0 + t1 * z) + lf(t0 ~ 1) + lf(t1 ~ 1), family = fm, data = dt)); cat("nlf(tau ~ t0 + t1 z)", ll(nl1), "\n")
nl2 <- tt(frm(bf(choice | reward(pay1, pay2) ~ a0 + a1 * z, a0 + a1 ~ 1, tau ~ 1, nl = TRUE), family = fm, data = dt)); cat("nl = TRUE (primary alpha)", ll(nl2), "\n")
l2 <- tt(frm(bf(choice | reward(pay1, pay2) ~ z, tau ~ 1), family = fm, data = dt)); cat("linear alpha ~ z", ll(l2), "\n")
nl3 <- tt(frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1) + nlf(alpha ~ a0 + a1 * z) + lf(a0 ~ 1) + lf(a1 ~ 1), family = fm, data = dt)); cat("nlf(alpha) ", ll(nl3), "\n")
if (!iserr(nl1)) {
  fr <- nl1$frame
  nll <- ns$build_objective(fr); tpl <- ns$make_start(fr, NULL)
  cat("learn guard elapsed:", system.time(ns$nl_check_body_read(fr, nll, tpl))[["elapsed"]], " fit elapsed:", system.time(suppressMessages(suppressWarnings(frm(bf(choice | reward(pay1, pay2) ~ 1) + nlf(tau ~ t0 + t1 * z) + lf(t0 ~ 1) + lf(t1 ~ 1), family = fm, data = dt))))[["elapsed"]], "\n")
}
cat("\n--- the table for the extension families\n")
for (f in c("cox", "royston_parmar", "wiener", "gddm", "lba", "rdm", "wiener_gng", "bandit2arm_delta", "rlddm", "hmm", "lca", "mixture_mvn")) {
  s <- tt(frm_compat("nl", f))
  if (iserr(s)) cat(f, "ERR", substr(s$msg, 1, 120), "\n") else { st <- as.data.frame(s); cat(sprintf("%-18s %s | %s\n", f, paste(st$status, collapse = ","), substr(paste(st$note %||% st$notes, collapse = " "), 1, 200))) }
}
cat("\n--- group:autocor x nl and the families the term fits under\n")
for (t in c("ar", "fcor", "cosy", "unstr")) { s <- tt(frm_compat(t, "nl")); if (!iserr(s)) { st <- as.data.frame(s); cat(t, paste(st$status, collapse = ","), "\n") } }
