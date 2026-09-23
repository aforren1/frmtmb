# Hands-on probe of drmTMB 0.7.0's model surface. Each probe prints the
# call and either a compact summary or the refusal message, so the log
# is the evidence and nothing here relies on drmTMB's own documentation.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(drmTMB))
cat("drmTMB", format(packageVersion("drmTMB")), "\n")

probe <- function(label, expr) {
  cl <- substitute(expr)
  cat("\n=== ", label, "\n", sep = "")
  cat(paste(deparse(cl, width.cutoff = 100), collapse = "\n"), "\n")
  res <- tryCatch(
    withCallingHandlers(eval(cl, parent.frame()),
                        warning = function(w) {
                          cat("  WARNING:", conditionMessage(w), "\n")
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) e)
  if (inherits(res, "error")) {
    msg <- gsub("\n", " ", conditionMessage(res))
    cat("  REFUSED/ERROR:", substr(msg, 1, 600), "\n")
    return(invisible(NULL))
  }
  if (inherits(res, "drmTMB")) {
    cat("  FITTED. logLik =", format(as.numeric(logLik(res)), digits = 12),
        " df =", attr(logLik(res), "df"), "\n")
    cf <- tryCatch(res$opt$par, error = function(e) NULL)
    if (!is.null(cf)) print(signif(cf, 6))
  } else {
    print(res)
  }
  invisible(res)
}

set.seed(20260922)
ng <- 30; nper <- 10; n <- ng * nper
g <- factor(rep(seq_len(ng), each = nper))
w <- rnorm(ng)[g]                    # group-level covariate
x <- rnorm(n); z <- rnorm(n)
u <- rnorm(ng, 0, 0.6)[g]; v <- rnorm(ng, 0, 0.3)[g]
d <- data.frame(g = g, w = w, x = x, z = z,
                y = 1 + 0.5 * x + u + rnorm(n, 0, exp(-0.2 + 0.3 * z + v)))
d$yb <- plogis(0.2 + 0.4 * x + u)
d$yb <- rbeta(n, d$yb * 20, (1 - d$yb) * 20)
d$yc <- rnbinom(n, mu = exp(1 + 0.3 * x + u), size = 3)
d$yp <- rpois(n, exp(0.5 + 0.3 * x + u))
d$ybin <- rbinom(n, 1, plogis(0.3 * x + u))
d$yt <- 1 + 0.5 * x + u + 0.8 * rt(n, df = 5)
d$yo <- cut(0.8 * x + u + rlogis(n), c(-Inf, -1, 0.5, 2, Inf),
            labels = FALSE)
d$yo <- factor(d$yo, ordered = TRUE)
d$y2 <- 0.5 * d$y + rnorm(n)
d$y3 <- rnorm(n)

cat("\n##### 1. Which distributional parameters take random effects\n")
probe("gaussian mu RE only",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d))
probe("gaussian sigma RE",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
             data = d))
probe("gaussian sigma random slope",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (0 + z | g)),
             family = gaussian(), data = d))
probe("gaussian correlated mu/sigma via label |p|",
      drmTMB(bf(y ~ x + (1 | p | g), sigma ~ z + (1 | p | g)),
             family = gaussian(), data = d))
probe("gaussian correlated intercept-slope block in mu",
      drmTMB(bf(y ~ x + (1 + x | g), sigma ~ 1), family = gaussian(),
             data = d))
probe("student mu RE",
      drmTMB(bf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = student(),
             data = d))
probe("student sigma RE",
      drmTMB(bf(yt ~ x, sigma ~ 1 + (1 | g), nu ~ 1), family = student(),
             data = d))
probe("student nu RE",
      drmTMB(bf(yt ~ x, sigma ~ 1, nu ~ 1 + (1 | g)), family = student(),
             data = d))
probe("beta mu RE",
      drmTMB(bf(yb ~ x + (1 | g), sigma ~ 1), family = beta(), data = d))
probe("beta sigma RE",
      drmTMB(bf(yb ~ x, sigma ~ 1 + (1 | g)), family = beta(), data = d))
probe("beta correlated mu slope block",
      drmTMB(bf(yb ~ x + (1 + x | g), sigma ~ 1), family = beta(), data = d))
probe("nbinom2 mu RE",
      drmTMB(bf(yc ~ x + (1 | g), sigma ~ 1), family = nbinom2(), data = d))
probe("nbinom2 sigma RE",
      drmTMB(bf(yc ~ x, sigma ~ 1 + (1 | g)), family = nbinom2(), data = d))
probe("nbinom2 mu RE + sigma RE together",
      drmTMB(bf(yc ~ x + (1 | g), sigma ~ 1 + (1 | g)), family = nbinom2(),
             data = d))
probe("cumulative_logit mu RE",
      drmTMB(bf(yo ~ x + (1 | g)), family = cumulative_logit(), data = d))
probe("poisson correlated slope block (1 + x | g)",
      drmTMB(bf(yp ~ x + (1 + x | g)), family = poisson(), data = d))

cat("\n##### 2. What sd(g) ~ x does\n")
probe("sd(g) ~ w, gaussian",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ w), family = gaussian(),
             data = d))
probe("sd(g) ~ x (row-level covariate, varies within group)",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ x), family = gaussian(),
             data = d))
probe("sd(g) ~ w, poisson",
      drmTMB(bf(yp ~ x + (1 | g), sd(g) ~ w), family = poisson(), data = d))
probe("sd(g) ~ w, nbinom2",
      drmTMB(bf(yc ~ x + (1 | g), sigma ~ 1, sd(g) ~ w), family = nbinom2(),
             data = d))
probe("sd(g) ~ w with REML",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ w), family = gaussian(),
             data = d, REML = TRUE))

cat("\n##### 3. REML scope\n")
probe("REML gaussian (1|g), sigma ~ 1",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ 1), family = gaussian(), data = d,
             REML = TRUE))
probe("REML gaussian (1|g), sigma ~ z",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
             REML = TRUE))
probe("REML gaussian sigma RE",
      drmTMB(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
             data = d, REML = TRUE))
probe("REML gaussian correlated (1 + x | g)",
      drmTMB(bf(y ~ x + (1 + x | g), sigma ~ 1), family = gaussian(),
             data = d, REML = TRUE))
probe("REML gaussian fixed only",
      drmTMB(bf(y ~ x, sigma ~ 1), family = gaussian(), data = d,
             REML = TRUE))
probe("REML student", drmTMB(bf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1),
                             family = student(), data = d, REML = TRUE))
probe("REML beta", drmTMB(bf(yb ~ x + (1 | g), sigma ~ 1), family = beta(),
                          data = d, REML = TRUE))
probe("REML nbinom2", drmTMB(bf(yc ~ x + (1 | g), sigma ~ 1),
                             family = nbinom2(), data = d, REML = TRUE))
probe("REML poisson", drmTMB(bf(yp ~ x + (1 | g)), family = poisson(),
                             data = d, REML = TRUE))
probe("REML binomial (1|g)", drmTMB(bf(ybin ~ x + (1 | g)),
                                    family = binomial(), data = d,
                                    REML = TRUE))
probe("REML cumulative_logit", drmTMB(bf(yo ~ x + (1 | g)),
                                      family = cumulative_logit(), data = d,
                                      REML = TRUE))
probe("REML biv_gaussian ordinary RE",
      drmTMB(bf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g)),
             family = biv_gaussian(), data = d, REML = TRUE))

cat("\n##### 4. Parameterizations and links (fixed-effect fits)\n")
show_pars <- function(fit) {
  print(signif(fit$opt$par, 7))
  cat("  fixef():\n"); print(tryCatch(fixef(fit), error = function(e) e))
}
f <- probe("nbinom2 sigma", drmTMB(bf(yc ~ x, sigma ~ 1), family = nbinom2(),
                                   data = d))
show_pars(f)
cat("  exp(-2 * log sigma) = theta:", exp(-2 * f$opt$par[["beta_sigma"]]),
    "  MASS::glm.nb theta:", MASS::glm.nb(yc ~ x, data = d)$theta, "\n")
f <- probe("beta sigma", drmTMB(bf(yb ~ x, sigma ~ 1), family = beta(),
                                data = d))
show_pars(f)
f <- probe("student nu link", drmTMB(bf(yt ~ x, sigma ~ 1, nu ~ 1),
                                     family = student(), data = d))
show_pars(f)
f <- probe("skew_normal", drmTMB(bf(yt ~ x, sigma ~ 1, nu ~ 1),
                                 family = skew_normal(), data = d))
show_pars(f)
f <- probe("cumulative_logit thresholds",
           drmTMB(bf(yo ~ x), family = cumulative_logit(), data = d))
show_pars(f)
print(summary(f))
f <- probe("biv_gaussian rho12 link",
           drmTMB(bf(mu1 = y ~ x, mu2 = y2 ~ x, rho12 = ~ 1),
                  family = biv_gaussian(), data = d))
show_pars(f)
print(summary(f))

cat("\n##### 5. Things frmtmb has: does drmTMB accept them?\n")
probe("nonlinear formula y ~ a * exp(b * x)",
      drmTMB(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
             family = gaussian(), data = d))
probe("nonlinear formula without nl flag",
      drmTMB(bf(y ~ exp(x)), family = gaussian(), data = d))
probe("three responses mvbind(y, y2, y3)",
      drmTMB(bf(mvbind(y, y2, y3) ~ x), family = gaussian(), data = d))
probe("three responses family c(gaussian x3)",
      drmTMB(bf(mu1 = y ~ x, mu2 = y2 ~ x, mu3 = y3 ~ x),
             family = c(gaussian(), gaussian(), gaussian()), data = d))
probe("mixed family bivariate gaussian + poisson",
      drmTMB(bf(mu1 = y ~ x, mu2 = yp ~ x),
             family = c(gaussian(), poisson()), data = d))
probe("mixed family bivariate gaussian + nbinom2",
      drmTMB(bf(mu1 = y ~ x, mu2 = yc ~ x),
             family = list(gaussian(), nbinom2()), data = d))
probe("|p| cross-formula correlation in nbinom2 mu and sigma",
      drmTMB(bf(yc ~ x + (1 | p | g), sigma ~ 1 + (1 | p | g)),
             family = nbinom2(), data = d))
probe("|p| cross-formula correlation in beta mu and sigma",
      drmTMB(bf(yb ~ x + (1 | p | g), sigma ~ 1 + (1 | p | g)),
             family = beta(), data = d))
probe("custom family object (quasi)",
      drmTMB(bf(y ~ x), family = quasi(), data = d))
probe("custom family list with loglik function",
      drmTMB(bf(y ~ x), family = structure(
        list(family = "mylaplace", link = "identity",
             loglik = function(y, mu, sigma) {
               -abs(y - mu) / sigma - log(2 * sigma)
             }),
        class = "family"), data = d))
probe("smooth s(x)", drmTMB(bf(y ~ s(x)), family = gaussian(), data = d))
probe("monotonic mo()", drmTMB(bf(y ~ mo(yo)), family = gaussian(), data = d))
probe("brms addition term y | trials()",
      drmTMB(bf(yc | trials(20) ~ x), family = binomial(), data = d))
probe("brms gr(g, cov = A) spelling",
      drmTMB(bf(y ~ x + (1 | gr(g, cov = A))), family = gaussian(), data = d))
probe("brms family object brms::student() stand-in: frmtmb-style family name",
      drmTMB(bf(y ~ x), family = "student", data = d))
probe("autocorrelation ar()",
      drmTMB(bf(y ~ x + ar(p = 1)), family = gaussian(), data = d))
probe("zero-inflation zi ~ on nbinom2", drmTMB(bf(yc ~ x, sigma ~ 1, zi ~ 1),
                                              family = nbinom2(), data = d))
probe("zero-inflation zi ~ (1|g)",
      drmTMB(bf(yc ~ x, sigma ~ 1, zi ~ 1 + (1 | g)), family = nbinom2(),
             data = d))

cat("\n##### 6. Methods and sampling\n")
f <- drmTMB(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d)
meth <- as.character(methods(class = "drmTMB"))
print(sort(unique(sub("[.]drmTMB$", "", meth))))
print(grep("sample|mcmc|stan|posterior|draw",
           ls(getNamespace("drmTMB")), value = TRUE, ignore.case = TRUE))
