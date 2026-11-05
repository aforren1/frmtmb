# Items 2, 3, 6: small measurements. Seeds per block.
source("dev/review-refuse-lib.R")
n <- 200; set.seed(10200); x <- rnorm(n)
# gradient: is 2.6e-4 and 1.5e-5 the optimizer's stop, or a density mismatch?
# Tighten nlminb and read the taped gradient; the L1 value identity is apart
set.seed(10202)
k <- rbinom(n, 1, 0.35)
y <- ifelse(k == 1, 2.5 + 0.5 * x, -2 - 0.3 * x) + rnorm(n)
ce <- ifelse(y > 3, 1L, ifelse(y < -3, -1L, 0L))
d <- data.frame(x = x, y = pmax(pmin(y, 3), -3), ce = ce)
tight <- frmtmb_control(optCtrl = list(iter.max = 5000, eval.max = 5000, rel.tol = 1e-14, x.tol = 1e-12))
go <- function(lab, form, fam, data) {
  a <- tr(frm(form, family = fam, data = data))
  b <- tr(frm(form, family = fam, data = data, control = tight))
  if (is_err(a) || is_err(b)) { pr("%-30s %s | %s", lab, msg(a, 90), msg(b, 90)); return(invisible()) }
  pr("%-30s default max|gr| %.2e, tightened %.2e; logLik %.8f %.8f", lab,
     max(abs(a$obj$gr(a$opt$par))), max(abs(b$obj$gr(b$opt$par))), as.numeric(logLik(a)), as.numeric(logLik(b)))
}
go("mixture, no cens", bf(y ~ 1, mu1 ~ x, mu2 ~ x), mixture(gaussian(), gaussian(), order = "none"), d)
go("mixture, cens", bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), mixture(gaussian(), gaussian(), order = "none"), d)
del <- 4 / sqrt(17)
ds <- data.frame(x = x, s = 0.1 + 0.4 * runif(n))
ds$y <- 0.5 + 0.7 * x + 0.8 * (del * abs(rnorm(n)) + sqrt(1 - del^2) * rnorm(n)) + rnorm(n, 0, ds$s)
go("skew_normal, no se", bf(y ~ x), skew_normal(), ds)
go("skew_normal, se", bf(y | se(s) ~ x), skew_normal(), ds)
# the skew_normal se() contribution by hand, brms's mean parameterization
ft <- frm(bf(y | se(s) ~ x), family = skew_normal(), data = ds)
print(unlist(ft$estimates[c("beta", "betad")]))
al <- ft$estimates$betad[["alpha_(Intercept)"]]
mu <- ft$estimates$beta[[1]] + ft$estimates$beta[[2]] * ds$x
dl <- al / sqrt(1 + al^2)
om <- ds$s / sqrt(1 - 2 * dl^2 / pi)
xi <- mu - om * dl * sqrt(2 / pi)
hand <- sum(log(2) - log(om) + dnorm((ds$y - xi) / om, log = TRUE) + pnorm(al * (ds$y - xi) / om, log.p = TRUE))
pr("skew_normal se(): logLik - hand (brms's mean parameterization) %.3e", as.numeric(logLik(ft)) - hand)

# the mixed-cor key through the survey's flow
set.seed(10203)
g <- factor(sample(1:15, n, TRUE))
d5 <- data.frame(x = x, g = g)
u <- rnorm(15); v <- 0.8 * u + rnorm(15, 0, 0.3)
d5$y <- 1 + 0.5 * x + 0.6 * u[g] + rnorm(n, 0, exp(-0.3 + 0.4 * v[g]))
for (sp in c("bf(y ~ x + (1 | gr(g, id = 'p')), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))",
             "bf(y ~ x + (1 | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p')))")) {
  r <- tr(meas(list(bf = sp, data = d5)))
  if (is_err(r)) { pr("%s ERROR %s", sp, msg(r)); next }
  pr("%-78s verdict=%s const=%s pts=[%s] grad=%s l1[%s]", sp, r$verdict, format(r$l1_const %||% NA, digits = 3),
     paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","), format(r$l1_grad %||% NA, digits = 3), substr(r$l1_error %||% "", 1, 100))
  a <- frm(eval(parse(text = sp), envir = asNamespace("frmtmb")), data = d5)
  b <- brms::make_standata(eval(parse(text = sp), envir = asNamespace("brms")), data = d5)
  co <- brms::make_stancode(eval(parse(text = sp), envir = asNamespace("brms")), data = d5)
  pr("   frm blocks %s, logLik %.4f; brms has a correlation: %s", paste(vapply(a$frame$re_blocks, function(k) paste0(k$covstruct, k$dim), ""), collapse = ","),
     as.numeric(logLik(a)), grepl("cor_1", co))
}
a <- frm(bf(y ~ x + (1 | gr(g, id = "p")), sigma ~ (1 | gr(g, id = "p"))), data = d5)
pr("with the correlation both ways: logLik %.4f", as.numeric(logLik(a)))

# mi() read by another response on the rows where that response is itself
# missing: unread per row; an all-missing mi() and rellib-r16's gaussian answer
d <- data.frame(x = x, y = rnorm(n), z = rnorm(n))
d$y <- NA_real_
a <- tr(frm(bf(y | mi() ~ x), data = d))
pr("gaussian, all missing, this build: %s", msg(a, 150))
b <- tr(brms::make_standata(brms::bf(y | mi() ~ x), data = d))
pr("  brms: %s", msg(b, 150))

# the positive custom family: valid_y spelled right, placeholder 0
fam_args <- list(family = "tfam", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dgamma(y, 2, dpars[["mu"]] / 2, log = TRUE)
  })
set.seed(10204)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.3 + 0.5 * x))); d$y[seq(5, 195, 10)] <- NA
print(names(formals(frmtmb_family)))
vy <- function(y, ...) if (any(y <= 0)) stop("tfam: response must be positive") else invisible(NULL)
fp <- do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi", valid_y = vy)))
r <- tr(frm(bf(y | mi() ~ x), family = fp, data = d))
pr("custom positive family, valid_y, mi, no chart: %s", msg(r, 180))
fp2 <- do.call(frmtmb_family, c(fam_args, list(accepts_aterms = "mi", valid_y = vy, mi_chart = "log")))
r <- tr(frm(bf(y | mi() ~ x), family = fp2, data = d))
pr("custom positive family, valid_y, mi, chart log: %s", msg(r, 180))
# a sibling that takes mi? any family of an extension
for (p in c("frmtmb.survival", "frmtmb.eam", "frmtmb.latent", "frmtmb.spatial", "frmtmb.learn", "frmtmb.ode", "frmtmb.coupling")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
fnm <- grep("mi_chart", unlist(lapply(c("frmtmb.survival", "frmtmb.eam", "frmtmb.latent"), function(p) ls(asNamespace(p)))), value = TRUE)
pr("extensions naming mi_chart: %d", length(fnm))
