# Lane optstart: the cold starting values (par_template()) of models
# where the least-squares start applies and where it does not, printed
# in hex so that two builds can be compared bit for bit
# (dev/optstart-guard-cmp.R). One process per library.
#   Rscript dev/optstart-guard.R <lib or ""> <out.tsv>
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(11)
n <- 60
d <- data.frame(x = rnorm(n), x2 = runif(n), g = factor(rep(1:10, 6)),
                w = runif(n, 0.5, 2), o = rnorm(n, 0, 0.1))
eta <- 1 + 3 * d$x
d$y <- 50 + 30 * d$x + rnorm(10, 0, 2)[d$g] + rnorm(n)
d$ypos <- exp(1 + 0.8 * d$x + rnorm(n, 0, 0.3))
d$ycnt <- rnbinom(n, mu = exp(1 + 0.8 * d$x), size = 3)
d$ypois <- rpois(n, exp(1 + 0.5 * d$x))
d$yprop <- plogis(0.5 + 1.2 * d$x + rnorm(n, 0, 0.3))
d$ybin <- rbinom(n, 1, plogis(d$x))
d$ycens <- ifelse(d$y > 60, 1, 0)
d$yord <- cut(d$y, c(-Inf, 30, 50, 70, Inf), ordered_result = TRUE)
d$yzero <- ifelse(runif(n) < 0.3, 0, d$y)
d$yhu <- ifelse(runif(n) < 0.3, 0, d$ypos)
d$se <- runif(n, 0.5, 1)
d$u <- runif(n, 1, 2)
for (k in 1:8) d[[paste0("p", k)]] <- d$u^k
d$ypoly <- sin(3 * d$u) + rnorm(n, 0, 0.2)
# an mi() response with a third of its values missing
d$ymi <- ifelse(runif(n) < 0.3, NA, d$y)
d$lo <- 0
d$hi <- 200
pt <- function(f, ...) {
  tryCatch(unclass(suppressMessages(par_template(f, data = d, ...))),
           error = function(e) paste("ERROR", conditionMessage(e)))
}
cases <- list(
  # where the rule does not apply: the start before this lane, bit for bit
  absent_poisson = pt(bf(ypois ~ x), family = poisson()),
  absent_binomial = pt(bf(ybin ~ x), family = bernoulli()),
  absent_icpt_only = pt(bf(y ~ 1 + (1 | g)), family = gaussian()),
  absent_nl = pt(bf(y ~ a + b * x, a ~ 1, b ~ 1, nl = TRUE),
                 family = gaussian()),
  absent_mixture = pt(bf(y ~ x), family = mixture(gaussian(), gaussian())),
  absent_cens = pt(bf(y | cens(ycens) ~ x), family = gaussian()),
  absent_trunc = pt(bf(y | trunc(lb = lo, ub = hi) ~ x), family = gaussian()),
  absent_se = pt(bf(y | se(se, sigma = TRUE) ~ x), family = gaussian()),
  absent_sigma_fixed = pt(bf(y ~ x, sigma = 2), family = gaussian()),
  absent_ordinal = pt(bf(yord ~ x), family = cumulative()),
  absent_gaussian_log = pt(bf(ypos ~ x), family = gaussian(link = "log")),
  absent_gamma_inverse = pt(bf(ypos ~ x), family = Gamma(link = "inverse")),
  absent_gamma_noicpt = pt(bf(ypos ~ 0 + x), family = Gamma()),
  absent_hurdle_gamma = pt(bf(yhu ~ x), family = hurdle_gamma()),
  absent_student_se = pt(bf(y | se(se, sigma = TRUE) ~ x), family = student()),
  absent_tweedie = pt(bf(ycnt ~ x), family = tweedie()),
  # a design whose condition number is above 1e5
  absent_illcond = pt(bf(ypoly ~ p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8),
                      family = gaussian()),
  # where it applies
  # a start or a prior: the rule applies to what start does not name
  start_theta = pt(bf(y ~ x + (1 | g)), family = gaussian(),
                   start = list(theta = 0)),
  prior_b = pt(bf(y ~ x + (1 | g)), family = gaussian(),
               prior = set_prior("normal(0, 100)", class = "b")),
  absent_smooth = pt(bf(y ~ s(x)), family = gaussian()),
  present_gaussian = pt(bf(y ~ x + (1 | g)), family = gaussian()),
  present_gaussian_w = pt(bf(y | weights(w) ~ x), family = gaussian()),
  present_gaussian_off = pt(bf(y ~ x + offset(o)), family = gaussian()),
  present_gaussian_sigx = pt(bf(y ~ x, sigma ~ x2), family = gaussian()),
  present_gaussian_noicpt = pt(bf(y ~ 0 + x), family = gaussian()),
  present_student = pt(bf(y ~ x), family = student()),
  present_skew_normal = pt(bf(y ~ x), family = skew_normal()),
  present_exgaussian = pt(bf(y ~ x), family = exgaussian()),
  present_asym_laplace = pt(bf(y ~ x), family = asym_laplace()),
  present_zi_asym_laplace = pt(bf(yzero ~ x),
                               family = zero_inflated_asym_laplace()),
  present_huber = pt(bf(y ~ x), family = huber()),
  present_lognormal = pt(bf(ypos ~ x), family = lognormal()),
  present_shifted_lognormal = pt(bf(ypos ~ x), family = shifted_lognormal()),
  present_hurdle_lognormal = pt(bf(yhu ~ x), family = hurdle_lognormal()),
  present_gamma = pt(bf(ypos ~ x), family = Gamma(link = "log")),
  present_weibull = pt(bf(ypos ~ x), family = weibull()),
  present_beta = pt(bf(yprop ~ x), family = Beta()),
  present_negbinomial = pt(bf(ycnt ~ x), family = negbinomial()),
  present_mi = pt(bf(ymi | mi() ~ x), family = gaussian()),
  present_mv = pt(mvbf(bf(y ~ x), bf(ypos ~ x, family = lognormal()),
                       rescor = FALSE), family = gaussian())
)
rows <- list()
for (nm in names(cases)) {
  v <- cases[[nm]]
  if (is.character(v)) {
    rows[[length(rows) + 1L]] <- data.frame(case = nm, par = "-", hex = v,
                                            value = NA)
    next
  }
  u <- unlist(v)
  rows[[length(rows) + 1L]] <- data.frame(
    case = nm, par = names(u), hex = sprintf("%a", u), value = u)
}
X <- do.call(rbind, rows)
write.table(X, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat("cases", length(cases), "values", nrow(X), "\n")
