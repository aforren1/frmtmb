# Item 6's bitwise battery. Every fit is recorded as logLik, opt$par and
# every standard error the fit reports (vcov(full = TRUE), vcov(), and
# the random effects' conditional sds); the arms are compared with
# identical() by dev/predfix-bitwise-cmp.R.
#   PREDFIX_ARM=base|lane PREDFIX_TAG=<tag> Rscript dev/predfix-bitwise.R
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
tag <- Sys.getenv("PREDFIX_TAG", PREDFIX_ARM)
set.seed(20260923)
n <- 200
x <- rnorm(n)
z <- runif(n)
g <- factor(rep(1:10, length.out = n))
f3 <- factor(sample(c("a", "b", "c"), n, TRUE))
u <- rnorm(10, 0, 0.4)[g]
eta <- 0.3 + 0.4 * x + u
pos <- function(v) exp(v)
d <- data.frame(x = x, z = z, g = g, f = f3, tr = 10L)
d$gau <- eta + rnorm(n)
d$pois <- rpois(n, exp(eta))
d$bern <- rbinom(n, 1, plogis(eta))
d$binom <- rbinom(n, 10, plogis(eta))
d$gam <- rgamma(n, 2, 2 / exp(eta))
d$lnorm <- exp(eta + rnorm(n, 0, 0.5))
d$stud <- eta + rt(n, 4)
d$nb <- rnbinom(n, mu = exp(eta), size = 3)
d$beta <- plogis(eta + rnorm(n, 0, 0.5))
d$zip <- ifelse(runif(n) < 0.2, 0L, rpois(n, exp(eta)))
d$zinb <- ifelse(runif(n) < 0.2, 0L, rnbinom(n, mu = exp(eta), size = 3))
d$hp <- ifelse(runif(n) < 0.3, 0L, 1L + rpois(n, exp(eta)))
d$geom <- rnbinom(n, mu = exp(eta), size = 1)
d$expo <- rexp(n, 1 / exp(eta))
d$weib <- rweibull(n, 2, exp(eta))
d$slnorm <- 0.2 + exp(eta + rnorm(n, 0, 0.5))
d$hgam <- ifelse(runif(n) < 0.3, 0, rgamma(n, 2, 2 / exp(eta)))
d$hlnorm <- ifelse(runif(n) < 0.3, 0, exp(eta + rnorm(n, 0, 0.5)))
d$zib <- ifelse(runif(n) < 0.2, 0L, rbinom(n, 10, plogis(eta)))
d$zibeta <- ifelse(runif(n) < 0.2, 0, plogis(eta + rnorm(n, 0, 0.5)))
d$alap <- eta + rnorm(n) * 0.7
d$hub <- eta + rt(n, 3)
d$bb <- rbinom(n, 10, plogis(eta + rnorm(n, 0, 0.5)))
d$skew <- eta + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
d$exg <- eta + rnorm(n, 0, 0.5) + rexp(n, 2)
d$ord <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
             labels = FALSE)
d$ordf <- factor(d$ord, ordered = TRUE)
d$cat <- factor(sample(c("p", "q", "r"), n, TRUE))
d$vm <- (0.3 * x + rnorm(n, 0, 0.5)) %% (2 * pi) - pi
d$tw <- ifelse(runif(n) < 0.2, 0, rgamma(n, 2, 2 / exp(eta)))

fams <- list(
  gaussian = list("gau", gaussian()), poisson = list("pois", poisson()),
  bernoulli = list("bern", bernoulli()),
  binomial = list("binom | trials(tr)", binomial()),
  Gamma = list("gam", Gamma(link = "log")),
  lognormal = list("lnorm", lognormal()), student = list("stud", student()),
  negbinomial = list("nb", negbinomial()), Beta = list("beta", Beta()),
  zero_inflated_poisson = list("zip", zero_inflated_poisson()),
  zero_inflated_negbinomial = list("zinb", zero_inflated_negbinomial()),
  hurdle_poisson = list("hp", hurdle_poisson()),
  geometric = list("geom", geometric()),
  exponential = list("expo", exponential()),
  weibull = list("weib", weibull()),
  shifted_lognormal = list("slnorm", shifted_lognormal()),
  hurdle_gamma = list("hgam", hurdle_gamma()),
  hurdle_lognormal = list("hlnorm", hurdle_lognormal()),
  zero_inflated_binomial = list("zib | trials(tr)", zero_inflated_binomial()),
  zero_inflated_beta = list("zibeta", zero_inflated_beta()),
  asym_laplace = list("alap", asym_laplace()), huber = list("hub", huber()),
  beta_binomial = list("bb | trials(tr)", beta_binomial()),
  skew_normal = list("skew", skew_normal()),
  exgaussian = list("exg", exgaussian()),
  cumulative = list("ordf", cumulative()), sratio = list("ord", sratio()),
  cratio = list("ord", cratio()), acat = list("ord", acat()),
  categorical = list("cat", categorical()),
  von_mises = list("vm", von_mises()), tweedie = list("tw", tweedie())
)
designs <- c(icpt = "1 + x", noint = "0 + x", noint_f = "0 + f + x",
             mixed = "x + (1 | g)", slope = "x + (1 + x | g)")
modes <- c("ML", "REML", "profile")

grab <- function(expr) tryCatch(expr, error = function(e) {
  structure(conditionMessage(e), class = "err")
})
record <- function(fit) {
  if (inherits(fit, "err")) return(list(error = unclass(fit)))
  list(ll = grab(as.numeric(logLik(fit))), par = fit$opt$par,
       conv = fit$opt$convergence,
       se_full = grab(sqrt(diag(vcov(fit, full = TRUE)))),
       se_fix = grab(sqrt(diag(vcov(fit)))),
       sd_re = grab(sqrt(frmtmb:::sdr_of(fit)$diag.cov.random)),
       units = fit$par_units)
}
fit_one <- function(fo, fam, data, mode, ctl = list()) {
  args <- list(bf(stats::as.formula(fo)), family = fam, data = data)
  if (mode == "REML") args$REML <- TRUE
  if (mode == "profile") ctl$profile <- TRUE
  args$control <- do.call(frmtmb_control, ctl)
  grab(suppressWarnings(suppressMessages(do.call(frm, args))))
}
out <- list()
t0 <- proc.time()[["elapsed"]]
for (fn in names(fams)) for (dn in names(designs)) for (m in modes) {
  y <- fams[[fn]][[1]]
  fo <- paste(y, "~", designs[[dn]])
  key <- paste(fn, dn, m, sep = "|")
  out[[key]] <- record(fit_one(fo, fams[[fn]][[2]], d, m))
}
# the targeted fits, where a column is spread below 1e-3 and the lane
# is MEANT to move (ML, or a dpar coefficient) or meant NOT to (a mu
# coefficient under REML or profile)
d$xs6 <- d$x * 1e-6
d$xs4 <- d$x * 1e-4
tgt <- list(
  "gaussian|1 + xs6" = list("gau ~ 1 + xs6", gaussian()),
  "gaussian|0 + xs6" = list("gau ~ 0 + xs6", gaussian()),
  "poisson|1 + xs6" = list("pois ~ 1 + xs6", poisson()),
  "poisson|0 + xs6" = list("pois ~ 0 + xs6", poisson()),
  "poisson|xs6 + (1 | g)" = list("pois ~ xs6 + (1 | g)", poisson()),
  "bernoulli|1 + xs6" = list("bern ~ 1 + xs6", bernoulli()),
  "negbinomial|1 + xs6" = list("nb ~ 1 + xs6", negbinomial()),
  "gaussian|1 + xs4" = list("gau ~ 1 + xs4", gaussian()),
  "poisson|0 + xs4" = list("pois ~ 0 + xs4", poisson()),
  "Gamma|1 + xs6" = list("gam ~ 1 + xs6", Gamma(link = "log"))
)
for (tn in names(tgt)) for (m in modes) {
  out[[paste("TARGET", tn, m, sep = "|")]] <-
    record(fit_one(tgt[[tn]][[1]], tgt[[tn]][[2]], d, m))
}
# a random slope on a column spread 1e-6 and on one spread 2e-2 (punch
# round 1, M1): engaged under every mode, because theta is outer
d$xs2 <- d$x * 2e-2
for (tn in c("xs6", "xs2")) for (m in modes) {
  out[[paste("TARGET", paste0("gaussian|", tn, " + (1 + ", tn, " | g)"), m,
             sep = "|")]] <- record(fit_one(
    paste0("gau ~ ", tn, " + (1 + ", tn, " | g)"), gaussian(), d, m))
}
# a dpar coefficient on a tiny column: outer under REML too
for (m in modes) {
  out[[paste("TARGET", "gaussian|x, sigma ~ xs6", m, sep = "|")]] <- record(
    grab(suppressWarnings(suppressMessages(frm(
      bf(gau ~ x, sigma ~ xs6), family = gaussian(), data = d,
      REML = m == "REML",
      control = frmtmb_control(profile = m == "profile"))))))
}
cat("fits:", length(out), " errors:",
    sum(vapply(out, function(r) !is.null(r$error), NA)),
    " seconds:", round(proc.time()[["elapsed"]] - t0), "\n")
saveRDS(out, paste0("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/",
                    "predfix-log/bitwise-", tag, ".rds"))
