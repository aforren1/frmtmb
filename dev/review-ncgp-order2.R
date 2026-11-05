# Item 4, fit side, second pass with the reader spellings this build takes: profile,
# bootstrap, hypothesis, allfit, sandwich, update, start, map, priors with coef; also
# readers after obj$fn moved the tape. Same metamorphic design as review-ncgp-order.R.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_tape_order
ident <- function(frame, tpl) tpl
setorder <- function(f) { unlockBinding("gp_tape_order", ns); assign("gp_tape_order", f, envir = ns); lockBinding("gp_tape_order", ns) }
set.seed(42)
n <- 60
d <- data.frame(x = runif(n, 0, 6), z = rnorm(n), g = factor(sample(letters[1:6], n, TRUE)), h = factor(sample(LETTERS[1:4], n, TRUE)),
                o = ordered(sample(1:4, n, TRUE)), f = factor(sample(c("u", "v"), n, TRUE)))
d$y <- 1 + sin(d$x) + 0.5 * d$z + rnorm(6, sd = .6)[as.integer(d$g)] + as.integer(d$o) * 0.3 + rnorm(n, sd = 0.25)
d$cnt <- rpois(n, exp(1 + 0.5 * sin(d$x) + 0.2 * d$z)); d$bin <- rbinom(n, 1, plogis(sin(d$x)))
d$y2 <- cos(d$x) + rnorm(n, sd = 0.3); d$zm <- d$z; d$zm[sample(n, 6)] <- NA
d$zobs <- d$z + rnorm(n, sd = 0.1); d$sdz <- 0.1
d$gv <- 2 + sin(d$x) + 0.4 * (((-log(runif(n)))^(-0.1) - 1) / 0.1)
models <- list(
  gp_re = list(f = y ~ z + gp(x) + (1 | g)),
  slope = list(f = y ~ z + gp(x) + (z | g) + (1 | h)),
  gp_s = list(f = y ~ gp(x) + s(z)),
  gp_mo = list(f = y ~ z + gp(x) + mo(o)),
  gp_me = list(f = y ~ me(zobs, sdz) + gp(x)),
  gp_mi = list(f = bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)),
  mvbf = list(f = bf(y ~ gp(x)) + bf(y2 ~ z + gp(x)) + set_rescor(FALSE)),
  mv_res = list(f = bf(y ~ gp(x)) + bf(y2 ~ z) + set_rescor(TRUE)),
  sigma_gp = list(f = bf(y ~ z, sigma ~ gp(x))),
  mix = list(f = bf(y ~ gp(x)), family = mixture(gaussian(), gaussian())),
  pois = list(f = cnt ~ z + gp(x) + (1 | g), family = poisson()),
  bern = list(f = bin ~ gp(x), family = bernoulli()),
  by = list(f = y ~ gp(x, by = f)),
  two = list(f = y ~ gp(x) + gp(z)),
  reml = list(f = y ~ z + gp(x), REML = TRUE),
  ordinal = list(f = o ~ gp(x), family = cumulative()),
  student = list(f = y ~ gp(x), family = student()),
  gev = list(f = gv ~ gp(x), family = gen_extreme_value()),
  map_pr = list(f = y ~ z + gp(x), prior = set_prior("normal(0, 1)", class = "b", coef = "z"))
)
pick <- if (length(args)) args else names(models)
tryr <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) structure(list(conditionMessage(err)), class = "err"))
readers <- function(fit, dat) {
  nm1 <- names(fixef(fit, flatten = TRUE))
  nm <- nm1[min(2L, length(nm1))]
  nd <- dat[1:7, , drop = FALSE]; nd$x <- nd$x + 0.0317
  list(
    ll = as.numeric(logLik(fit)),
    fixef = tryr(fixef(fit, flatten = TRUE)), coef = tryr(unlist(coef(fit))), ranef = tryr(unlist(ranef(fit))),
    vcov = tryr(as.matrix(vcov(fit))), vcovn = tryr(rownames(as.matrix(vcov(fit)))), vcovf = tryr(as.matrix(vcov(fit, full = TRUE))),
    ci = tryr(confint(fit)), prof = tryr(confint(fit, method = "profile", parm = nm)),
    sum = tryr(summary(fit)$coefficients), vc = tryr(unlist(VarCorr(fit))),
    fitted = tryr(fitted(fit)), fittedn = tryr(fitted(fit, newdata = nd)), linp = tryr(frm_linpred(fit, newdata = nd, se.fit = TRUE)),
    resid = tryr(residuals(fit)), sim = tryr({ set.seed(9); unlist(simulate(fit, nsim = 2, seed = 1)) }),
    getb = tryr(frmtmb:::getME.frmtmb_fit(fit, "b")), hc = tryr(as.matrix(vcov(fit, type = "HC0"))), cr = tryr(as.matrix(vcov(fit, cluster = dat$g))),
    start_th = tryr(as.numeric(logLik(update(fit, start = list(theta = fit$estimates$theta))))),
    start_b = tryr(as.numeric(logLik(update(fit, start = list(b = fit$estimates$b))))),
    upd = tryr(as.numeric(logLik(update(fit)))),
    boot = tryr({ b <- frm_bootstrap(fit, nsim = 2, seed = 3); b$t }),
    allfit = tryr({ a <- frm_allfit(fit, optimizers = "nlminb"); as.numeric(unlist(a$llik %||% a$summary$logLik %||% a[[1]])) }),
    hyp = tryr(hypothesis(fit, paste0(nm, " > 0"))$hypothesis[, c("Estimate", "Est.Error")]),
    diag = tryr(unlist(diagnose(fit)$summary)),
    parT = tryr(unlist(par_template(fit))),
    mapfit = tryr({ m2 <- update(fit, map = list(theta = factor(c(NA, seq_along(fit$estimates$theta)[-1])))); as.numeric(logLik(m2)) })
  )
}
num <- function(x) suppressWarnings(as.numeric(unlist(x)))
cmp <- function(a, b) {
  if (inherits(a, "err") || inherits(b, "err")) {
    if (inherits(a, "err") && inherits(b, "err")) return(if (identical(a[[1]], b[[1]])) paste0("both-refuse:", substr(a[[1]], 1, 60)) else "DIFF-ERR")
    return(paste0("ONE-ERR:", substr(if (inherits(a, "err")) a[[1]] else b[[1]], 1, 110)))
  }
  if (is.character(a) || is.character(b)) return(if (identical(a, b)) "same-names" else if (setequal(a, b)) "PERM-names" else "DIFF-names")
  na <- num(a); nb <- num(b)
  if (length(na) != length(nb)) return(sprintf("LEN %d/%d", length(na), length(nb)))
  if (!length(na)) return("empty")
  if (any(is.na(na) != is.na(nb))) return(sprintf("NA-PATTERN %d/%d", sum(is.na(na)), sum(is.na(nb))))
  ok <- !is.na(na)
  if (!any(ok)) return("allNA")
  sprintf("%.1e", max(abs(na[ok] - nb[ok]) / pmax(abs(nb[ok]), 1e-3)))
}
for (mn in pick) {
  m <- models[[mn]]; a <- m; a$f <- NULL
  setorder(orig)
  fa <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(m$f, data = d), a)))), error = function(e) e)
  setorder(ident)
  fb <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(m$f, data = d), a)))), error = function(e) e)
  if (inherits(fa, "error") || inherits(fb, "error")) {
    cat("RO", mn, "FIT ERROR lane:", if (inherits(fa, "error")) substr(conditionMessage(fa), 1, 140) else "ok", "| ident:", if (inherits(fb, "error")) substr(conditionMessage(fb), 1, 140) else "ok", "\n"); next
  }
  setorder(orig); ra <- readers(fa, d)
  setorder(ident); rb <- readers(fb, d)
  setorder(orig)
  out <- vapply(names(ra), function(k) cmp(ra[[k]], rb[[k]]), "")
  cat(sprintf("RO %s order %s | ident %s | conv %d/%d\n", mn, paste(names(fa$obj$env$parameters), collapse = ","), paste(names(fb$obj$env$parameters), collapse = ","), fa$opt$convergence, fb$opt$convergence))
  for (k in names(out)) cat(sprintf("RO   %-8s %s\n", k, out[[k]]))
  # a moved tape: what fitted() and summary() read after the tape moved
  se0 <- tryr(fitted(fa)[, 2]); v0 <- tryr(as.matrix(vcov(fa)))
  q <- fa$opt$par + 0.2
  fa$obj$fn(q); fa$obj$gr(q)
  se1 <- tryr(fitted(fa)[, 2]); v1 <- tryr(as.matrix(vcov(fa)))
  cat(sprintf("RO   moved tape: fitted se %s vcov %s\n", cmp(se1, se0), cmp(v1, v0)))
  flush.console()
}
cat("DONE\n")
