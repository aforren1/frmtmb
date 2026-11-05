# Item 4, fit side: the metamorphic test. Fit each model with b first (the lane)
# and with the tape in template order (gp_tape_order() patched to the identity in
# this process), and read every positional reader on both. Any difference above
# rounding is a reader keyed on position.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_tape_order
ident <- function(frame, tpl) tpl
setorder <- function(f) {
  unlockBinding("gp_tape_order", ns)
  assign("gp_tape_order", f, envir = ns)
  lockBinding("gp_tape_order", ns)
}
set.seed(42)
n <- 90
d <- data.frame(x = runif(n, 0, 6), z = rnorm(n), g = factor(sample(letters[1:8], n, TRUE)),
                h = factor(sample(LETTERS[1:4], n, TRUE)), o = ordered(sample(1:4, n, TRUE)),
                f = factor(sample(c("u", "v"), n, TRUE)), w = runif(n, 0.5, 2))
re <- rnorm(8, sd = 0.6)
d$y <- sin(d$x) + 0.5 * d$z + re[as.integer(d$g)] + as.integer(d$o) * 0.3 + rnorm(n, sd = 0.25)
d$cnt <- rpois(n, exp(1 + 0.5 * sin(d$x) + 0.2 * d$z))
d$bin <- rbinom(n, 1, plogis(sin(d$x)))
d$y2 <- cos(d$x) + rnorm(n, sd = 0.3)
d$zm <- d$z; d$zm[sample(n, 10)] <- NA
d$zobs <- d$z + rnorm(n, sd = 0.1)
d$sdz <- 0.1
models <- list(
  gp_re = list(f = y ~ z + gp(x) + (1 | g)),
  re_slope = list(f = y ~ z + gp(x) + (z | g) + (1 | h)),
  gp_s = list(f = y ~ gp(x) + s(z)),
  gp_mo = list(f = y ~ z + gp(x) + mo(o)),
  gp_me = list(f = y ~ me(zobs, sdz) + gp(x)),
  gp_mi = list(f = bf(y ~ mi(zm) + gp(x)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)),
  mvbf = list(f = bf(y ~ gp(x)) + bf(y2 ~ z + gp(x)) + set_rescor(FALSE)),
  mvbf_res = list(f = bf(y ~ gp(x)) + bf(y2 ~ z) + set_rescor(TRUE)),
  sigma_gp = list(f = bf(y ~ z, sigma ~ gp(x))),
  mix = list(f = bf(y ~ gp(x)), family = mixture(gaussian(), gaussian())),
  pois = list(f = cnt ~ z + gp(x) + (1 | g), family = poisson()),
  bern = list(f = bin ~ gp(x), family = bernoulli()),
  by = list(f = y ~ gp(x, by = f)),
  two = list(f = y ~ gp(x) + gp(z)),
  reml = list(f = y ~ z + gp(x), REML = TRUE),
  ordinal = list(f = o ~ gp(x), family = cumulative()),
  student = list(f = y ~ gp(x), family = student()),
  ar = list(f = y ~ gp(x) + ar(gr = g), start = NULL)
)
which_models <- if (length(args)) args else names(models)
num <- function(x) as.numeric(unlist(x))
tryr <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(err) structure(list(conditionMessage(err)), class = "err"))
readers <- function(fit, dat) {
  nd <- dat[1:7, , drop = FALSE]
  nd$x <- nd$x + 0.0317
  list(
    logLik = as.numeric(logLik(fit)),
    coefp = tryr(names(fit$opt$par)),
    fixef = tryr(fixef(fit)),
    coef = tryr(coef(fit)),
    ranef = tryr(ranef(fit)),
    vcov = tryr(vcov(fit)),
    vcovn = tryr(rownames(as.matrix(vcov(fit)))),
    confint = tryr(confint(fit)),
    sum = tryr(summary(fit)$coefficients),
    vc = tryr(VarCorr(fit)),
    fitted = tryr(fitted(fit)),
    pred = tryr(predict(fit, newdata = nd)),
    fittedn = tryr(fitted(fit, newdata = nd)),
    resid = tryr(residuals(fit)),
    prof = tryr(confint(fit, method = "profile", parm = names(fit$opt$par)[1])),
    gmeb = tryr(getME(fit, "b")),
    sim = tryr({ set.seed(9); simulate(fit, nsim = 2) }),
    sw = tryr(vcov(fit, type = "HC0")),
    start = tryr({ fit2 <- update(fit, start = list(theta = fit$estimates$theta), recompile = FALSE); as.numeric(logLik(fit2)) }),
    upd = tryr({ fit3 <- update(fit); as.numeric(logLik(fit3)) }),
    boot = tryr({ set.seed(3); b <- frm_bootstrap(fit, nsim = 2); as.numeric(unlist(b$estimates %||% b$t %||% b[[1]])) }),
    allfit = tryr(as.numeric(unlist(lapply(frm_allfit(fit, optimizers = "nlminb"), function(z) as.numeric(logLik(z)))))),
    hyp = tryr(hypothesis(fit, paste0(names(fixef(fit))[1], " > 0"))$hypothesis$Estimate),
    joint = tryr(dim(frmtmb:::get_joint_cov(fit)$V %||% matrix(0))),
    bse = tryr(frmtmb:::get_joint_cov(fit)$names),
    autosc = tryr(frmtmb:::autoscale_sdreport(fit)$cov.fixed)
  )
}
cmp <- function(a, b) {
  if (inherits(a, "err") || inherits(b, "err")) {
    if (inherits(a, "err") && inherits(b, "err")) return(if (identical(a[[1]], b[[1]])) "both-refuse" else "DIFF-ERR")
    return(paste0("ONE-ERR:", substr(if (inherits(a, "err")) a[[1]] else b[[1]], 1, 90)))
  }
  if (is.character(a) || is.character(b)) return(if (identical(a, b)) "same-names" else if (setequal(a, b)) "PERM-names" else "DIFF-names")
  na <- tryCatch(num(a), error = function(e) NULL); nb <- tryCatch(num(b), error = function(e) NULL)
  if (is.null(na) || is.null(nb)) return("?")
  if (length(na) != length(nb)) return(sprintf("LEN %d/%d", length(na), length(nb)))
  if (!length(na)) return("empty")
  if (all(is.na(na) == is.na(nb)) ) {
    ok <- !is.na(na)
    if (!any(ok)) return("allNA")
    return(sprintf("%.1e", max(abs(na[ok] - nb[ok]) / pmax(abs(nb[ok]), 1e-6))))
  }
  "NA-PATTERN"
}
for (mn in which_models) {
  m <- models[[mn]]
  a <- m; a$f <- NULL; a$start <- NULL
  setorder(orig)
  fa <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(m$f, data = d), a)))), error = function(e) e)
  setorder(ident)
  fb <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(list(m$f, data = d), a)))), error = function(e) e)
  if (inherits(fa, "error") || inherits(fb, "error")) {
    cat("RO", mn, "FIT ERROR lane:", if (inherits(fa, "error")) substr(conditionMessage(fa), 1, 140) else "ok",
        "| ident:", if (inherits(fb, "error")) substr(conditionMessage(fb), 1, 140) else "ok", "\n")
    next
  }
  setorder(orig); ra <- readers(fa, d)
  setorder(ident); rb <- readers(fb, d)
  setorder(orig)
  out <- vapply(names(ra), function(k) cmp(ra[[k]], rb[[k]]), "")
  cat(sprintf("RO %s tape-order %s | ident %s | conv %d/%d\n", mn, paste(names(fa$obj$env$parameters), collapse = ","),
              paste(names(fb$obj$env$parameters), collapse = ","), fa$opt$convergence, fb$opt$convergence))
  for (k in names(out)) cat(sprintf("RO   %-8s %s\n", k, out[[k]]))
  flush.console()
}
cat("DONE\n")
