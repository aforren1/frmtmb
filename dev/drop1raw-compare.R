# frmtmb's drop1() against stats::drop1() on lm() and glm() and lme4's
# drop1() on lmer(REML = FALSE), on models with transformed terms.
# Arg: base | lane. Log: dev/drop1raw-out/compare-<arm>.txt
# Data: seed 5, n = 200 (the survey's), plus three NA cells in x3.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
cat("lme4", format(packageVersion("lme4")), " brms",
    format(packageVersion("brms")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 200
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + 0.3 * d$x * log(d$x2) + rnorm(8)[d$g] +
  rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x + 0.3 * log(d$x2)))
dn <- d
dn$x3[c(3, 17, 40)] <- NA

# frmtmb's table: LRT and Df per dropped term
ft <- function(f, data, fam = gaussian(), scope = NULL) {
  fit <- q(frm(bf(f), family = fam, data = data))
  t <- if (is.null(scope)) q(drop1(fit, test = "Chisq")) else
    q(drop1(fit, scope = scope, test = "Chisq"))
  list(rows = rownames(t)[-1L], lrt = t$LRT[-1L], df = t$Df[-1L])
}
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1))
line <- function(lab, a, r) {
  if (inherits(a, "error")) {
    cat(sprintf("%-34s ERROR %s\n", lab, substr(conditionMessage(a), 1, 70)))
    return(invisible())
  }
  if (inherits(r, "error")) {
    cat(sprintf("%-34s REF ERROR %s\n", lab,
                substr(conditionMessage(r), 1, 70)))
    return(invisible())
  }
  cat(sprintf("%-34s rows %s | df %s | LRT max rel diff %.2e | LRT %s\n",
              lab, identical(a$rows, r$rows),
              identical(as.numeric(a$df), as.numeric(r$df)),
              rel(a$lrt, r$lrt),
              paste(sprintf("%.6f", a$lrt), collapse = " ")))
}

cat("\n== gaussian against lm(); LRT = n log(RSS_i / RSS_0) ==\n")
ref_lm <- function(f, data, scope = NULL) {
  m <- lm(f, data = data)
  t <- if (is.null(scope)) drop1(m) else drop1(m, scope = scope)
  nn <- nobs(m)
  list(rows = rownames(t)[-1L], df = t$Df[-1L],
       lrt = nn * log(t$RSS[-1L] / t$RSS[1L]))
}
lmcases <- list(
  log = z ~ x + log(x2), scale = z ~ scale(x) + x2, I = z ~ I(x^2) + x2,
  ns = z ~ x + splines::ns(x2, 3), poly_x = z ~ poly(x2, 2) + x,
  inter = z ~ x * log(x2), resp = log(z + 10) ~ x + log(x2))
for (nm in names(lmcases)) {
  f <- lmcases[[nm]]
  line(paste("lm", nm), tryCatch(ft(f, d), error = function(e) e),
       tryCatch(ref_lm(f, d), error = function(e) e))
}
line("lm scope = ~ x", tryCatch(ft(z ~ x + log(x2) + x3, d, scope = ~ x),
                                error = function(e) e),
     tryCatch(ref_lm(z ~ x + log(x2) + x3, d, scope = ~ x),
              error = function(e) e))
line("lm NA in x3, ns()", tryCatch(ft(z ~ splines::ns(x2, 3) + x3, dn),
                                   error = function(e) e),
     tryCatch(ref_lm(z ~ splines::ns(x2, 3) + x3, dn),
              error = function(e) e))
line("lm NA in x3, scale()", tryCatch(ft(z ~ scale(x) + x3, dn),
                                      error = function(e) e),
     tryCatch(ref_lm(z ~ scale(x) + x3, dn), error = function(e) e))
line("lm NA in x3, poly()", tryCatch(ft(z ~ poly(x2, 2) + x3, dn),
                                     error = function(e) e),
     tryCatch(ref_lm(z ~ poly(x2, 2) + x3, dn), error = function(e) e))

cat("\n== poisson against glm(); LRT = deviance difference ==\n")
ref_glm <- function(f, data) {
  m <- glm(f, family = poisson(), data = data)
  t <- drop1(m, test = "Chisq")
  list(rows = rownames(t)[-1L], df = t$Df[-1L], lrt = t$LRT[-1L])
}
glmcases <- list(log = cnt ~ x + log(x2),
                 scale_ns = cnt ~ scale(x) + splines::ns(x2, 3),
                 inter = cnt ~ x * log(x2))
for (nm in names(glmcases)) {
  f <- glmcases[[nm]]
  line(paste("glm", nm),
       tryCatch(ft(f, d, poisson()), error = function(e) e),
       tryCatch(ref_glm(f, d), error = function(e) e))
}
line("glm NA in x3, ns()",
     tryCatch(ft(cnt ~ splines::ns(x2, 3) + x3, dn, poisson()),
              error = function(e) e),
     tryCatch(ref_glm(cnt ~ splines::ns(x2, 3) + x3, dn),
              error = function(e) e))

cat("\n== gaussian mixed against lme4 lmer(REML = FALSE) ==\n")
# lme4's drop1() re-evaluates the lmer() call in the formula's
# environment, so the call has to name the data by a global symbol: a
# function argument `data` resolves there to utils::data() and the
# refits fail with "first argument must be a named list"
ref_lmer <- function(f, data) {
  dname <- deparse(substitute(data))
  m <- eval(bquote(lme4::lmer(.(f), data = .(as.name(dname)),
                              REML = FALSE)), globalenv())
  t <- drop1(m, test = "Chisq")
  list(rows = rownames(t)[-1L], df = t$npar[-1L], lrt = t$LRT[-1L])
}
lmercases <- list(log = z ~ x + log(x2) + (1 | g),
                  scale_ns = z ~ scale(x) + splines::ns(x2, 3) + (1 | g),
                  inter = z ~ x * log(x2) + (1 | g),
                  poly_x = z ~ poly(x2, 2) + x + (1 | g))
for (nm in names(lmercases)) {
  f <- lmercases[[nm]]
  line(paste("lmer", nm), tryCatch(ft(f, d), error = function(e) e),
       tryCatch(ref_lmer(f, d), error = function(e) e))
}
line("lmer NA in x3, ns()",
     tryCatch(ft(z ~ splines::ns(x2, 3) + x3 + (1 | g), dn),
              error = function(e) e),
     tryCatch(ref_lmer(z ~ splines::ns(x2, 3) + x3 + (1 | g), dn),
              error = function(e) e))

cat("\n== brms ==\n")
ns_brms <- asNamespace("brms")
cat("brms has a drop1 method:",
    any(grepl("^drop1[.]", ls(ns_brms, all.names = TRUE))) ||
      !is.null(utils::getS3method("drop1", "brmsfit", optional = TRUE)),
    "\n")
