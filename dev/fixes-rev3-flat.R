# Reviewer of lane fixes, final check: nl_flat_message() on models the
# lane did not run. The step is 1e-4 * max(|p|, 1) in raw coefficient
# units, so a coefficient of a covariate in large units moves the body a
# long way; many nonlinear coefficients and random effects in the
# nlpars; an optimum at a bound; and an unidentified ridge whose
# differenced Hessian is non-finite. Each line: the flat warning (if
# any), the smallest eigenvalue ratio the check computes, and the fit.
#   Rscript dev/fixes-rev3-flat.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
has_check <- exists("nl_flat_message", ns)
orig <- if (has_check) get("nl_flat_message", ns) else function() { flat <- ratio < nl_flat_tol }
txt <- deparse(orig)
hit <- grep("flat <- ratio < nl_flat_tol", txt, fixed = TRUE)
stopifnot(length(hit) == 1L)
txt[hit] <- "return(min(ratio))"
ratio_of <- eval(parse(text = txt))
environment(ratio_of) <- ns
run <- function(lab, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers({
    fit <- expr
    fe <- fixef(fit)
    rr <- tryCatch(ratio_of(fit$obj, fit$opt, fit$frame), error = function(e) NA)
    sprintf("ratio %s conv=%d est %s se %s", if (is.null(rr)) "NULL (H non-finite or < 2 coefs)" else format(rr, digits = 3),
            fit$opt$convergence,
            paste(sprintf("%.4g", fe[, "Estimate"]), collapse = " "),
            paste(sprintf("%.3g", fe[, "Est.Error"]), collapse = " "))
  }, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 150)))
  fl <- grep("are not identified", w, value = TRUE)
  cat(sprintf("%-44s FLAT %s | %s\n   other warnings: %s\n", lab,
              if (length(fl)) substr(sub(" of the nonlinear.*", "", fl[1]),
                                     1, 90) else "none", r,
              if (length(setdiff(w, fl))) substr(paste(setdiff(w, fl),
                                                       collapse = " / "),
                                                 1, 160) else "none"))
}
# 1. a covariate in large units: b ~ 2e-5, the step moves b * x by 10
set.seed(41)
x <- runif(300, 0, 1e5)
d1 <- data.frame(x = x, y = 3 * exp(-2e-5 * x) + rnorm(300, 0, 0.05))
run("a * exp(b * x), x 0..1e5, b ~ -2e-5",
    frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d1,
        start = list(beta = c(3, -2e-5))))
run("same, unidentified a1 + a2 shared",
    frm(bf(y ~ (a1 + a2) * exp(b * x), a1 ~ 1, a2 ~ 1, b ~ 1, nl = TRUE),
        data = d1, start = list(beta = c(1.5, 1.5, -2e-5))))
# 2. many nonlinear coefficients: a ~ 0 + f (40 levels), b ~ 1 + 3 covs
set.seed(42)
n <- 2000
d2 <- data.frame(f = factor(sample(1:40, n, TRUE)), z1 = rnorm(n),
                 z2 = rnorm(n), z3 = rnorm(n), x = runif(n))
af <- rnorm(40, 2, 0.4)
d2$y <- af[d2$f] * exp((0.5 + 0.2 * d2$z1 - 0.1 * d2$z2 + 0.05 * d2$z3) *
                         d2$x) + rnorm(n, 0, 0.3)
run("a ~ 0 + f (40), b ~ 1 + z1 + z2 + z3",
    frm(bf(y ~ a * exp(b * x), a ~ 0 + f, b ~ 1 + z1 + z2 + z3, nl = TRUE),
        data = d2, start = list(beta = c(rep(2, 40), 0.5, 0, 0, 0))))
# near-collinear covariates in b (cor 0.999)
d2$z4 <- d2$z1 + rnorm(n, 0, 0.045)
run("b ~ 1 + z1 + z4 (cor 0.999), a ~ 0 + f",
    frm(bf(y ~ a * exp(b * x), a ~ 0 + f, b ~ 1 + z1 + z4, nl = TRUE),
        data = d2, start = list(beta = c(rep(2, 40), 0.5, 0.1, 0.1))))
# 3. random effects in both nlpars, 30 groups
set.seed(43)
d3 <- data.frame(g = factor(rep(1:30, each = 20)), x = runif(600, 0, 3))
ra <- rnorm(30, 0, 0.3); rb <- rnorm(30, 0, 0.1)
d3$y <- (2 + ra[d3$g]) * exp((0.4 + rb[d3$g]) * d3$x) + rnorm(600, 0, 0.3)
run("a ~ 1 + (1 | g), b ~ 1 + (1 | g)",
    frm(bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1 + (1 | g), nl = TRUE),
        data = d3, start = list(beta = c(2, 0.4))))
# 4. an optimum at a bound (prior bound lb = 0 on a coefficient whose
# ML value is negative)
set.seed(44)
d4 <- data.frame(x = runif(200, 0, 3))
d4$y <- 2 * exp(-0.3 * d4$x) + rnorm(200, 0, 0.2)
run("b at its bound lb = 0",
    frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d4,
        start = list(beta = c(2, 0.1)),
        prior = set_prior("normal(0, 5)", nlpar = "b", lb = 0)))
# 5. an exact linear ridge whose differenced Hessian is non-finite: c0
# sits near 0, so c0 - h < 0 and log() is NaN there
set.seed(45)
d5 <- data.frame(x = rnorm(200))
d5$y <- 1 + 0.5 * d5$x + log(5e-5) + rnorm(200, 0, 0.3)
run("a + b + log(c0), a, b ~ 1 + x, c0 ~ 1 (c0 near 0)",
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = d5, start = list(beta = c(0.5, 0.25, 0.5, 0.25, 5e-5))))
run("same ridge, c0 near 1 (control)",
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = transform(d5, y = y - log(5e-5)),
        start = list(beta = c(0.5, 0.25, 0.5, 0.25, 1))))
# 6. update() and refit() of an unidentified fit: does the user see it?
d6 <- data.frame(x = rnorm(100)); d6$y <- 1 + d6$x + rnorm(100)
f6 <- suppressWarnings(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 0 + x, nl = TRUE),
                           data = d6))
run("update() to the ridge b ~ 1 + x",
    update(f6, formula. = bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE)))
f6b <- suppressWarnings(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                            data = d6))
run("refit() of the ridge to a new response",
    refit(f6b, d6$y + rnorm(100, 0, 0.1)))
