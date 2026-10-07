# Reviewer of lane setier: copy of dev/nanse-rev-cases.R with lane = wt-setier-lib over
# rellib-r6, base = rellib-r6 (merge arm = rellib-r6 alone, i.e. base).
# Reviewer: adversarial cases for the standard-error repair of lane
# nanse. Per case: optimizer code, every warning raised by frm() and by
# vcov(), the reported SEs, the lost set, and a reference.
#   Rscript dev/nanse-rev-cases.R base|lane|merge [case ...]
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  merge = c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
want <- if (length(args) > 1) args[-1] else NULL
on <- function(k) is.null(want) || k %in% want
ns <- asNamespace("frmtmb")

run <- function(lab, expr, ref = NULL) {
  cat("\n==", lab, "\n")
  w <- character()
  f <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  wv <- character()
  sdr <- withCallingHandlers(ns$sdr_of(f), warning = function(x) {
    wv <<- c(wv, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  V <- sdr$cov.fixed
  se <- suppressWarnings(sqrt(diag(V)))
  names(se) <- ns$outer_par_names(f)
  cat("  code", f$opt$convergence, "logLik", format(logLik(f), digits = 10),
      "\n")
  cat("  frm() warnings:", if (length(w)) {
    paste0("\n    ", substr(w, 1, 400), collapse = "")
  } else "none", "\n")
  cat("  sdr_of() warnings:", if (length(wv)) {
    paste0("\n    ", substr(wv, 1, 200), collapse = "")
  } else "none", "\n")
  show <- if (length(se) > 12) {
    c(head(se, 6), tail(se, 6))
  } else se
  cat("  SE:", paste(sprintf("%s=%.4g", names(show), show), collapse = " "),
      if (length(se) > 12) paste0(" [", length(se), " in all]"), "\n")
  cat("  lost:", if (length(sdr$se_lost)) {
    paste(sprintf("%s(%s)", names(sdr$se_lost), sdr$se_lost),
          collapse = " ")
  } else "none", "\n")
  if (!is.null(ref)) {
    r <- ref(f, se)
  }
  invisible(list(fit = f, se = se, warn = w))
}

# ---- 1. a null direction spread over many parameters ------------------
# y ~ a + b with a ~ 0 + f: only a_k + b is identified. In unit-diagonal
# coordinates the null vector loads sqrt(1/(2k)) on each a_k and 0.707
# on b, so past k = 50 no a_k reaches se_load_tol = 0.1.
spread <- function(k, m = 5, seed = 1) {
  set.seed(seed)
  f <- factor(rep(seq_len(k), each = m))
  mu <- rnorm(k)[f]
  data.frame(f = f, y = mu + rnorm(k * m, 0, 0.5))
}
ref_ridge <- function(f, se) {
  # the identified quantity is a_k + b; the SE of a_k alone is infinite
  nm <- names(se)
  cat("  finite SEs on a_*:", sum(is.finite(se[grepl("^a_", nm)])), "of",
      sum(grepl("^a_", nm)), "; on b_*:", sum(is.finite(se[grepl("^b_", nm)])),
      "\n")
}
for (k in c(10L, 40L, 60L, 120L)) {
  if (!on(paste0("spread", k))) next
  dd <- spread(k)
  run(paste0("spread ridge y ~ a + b, a ~ 0 + f, k = ", k),
      frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE), data = dd,
          family = gaussian()), ref = ref_ridge)
}

# ---- 2. perfectly separated logistic regression ----------------------
if (on("sep")) {
  set.seed(2)
  n <- 60
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- as.integer(dd$x > 0)
  g <- suppressWarnings(glm(y ~ x + z, binomial, data = dd))
  cat("\n glm: coef", paste(signif(coef(g), 4), collapse = " "), " SE",
      paste(signif(sqrt(diag(vcov(g))), 4), collapse = " "), "\n")
  run("separated logistic y ~ x + z, y = x > 0",
      frm(bf(y ~ x + z), family = bernoulli(), data = dd))
}

# ---- 3. a variance component at 0 --------------------------------------
if (on("vc0")) {
  cnt <- c(tier = 0, warn = 0, nan = 0)
  for (s in 1:20) {
    set.seed(100 + s)
    dd <- data.frame(g = factor(rep(1:8, each = 6)), x = rnorm(48))
    dd$y <- 1 + 0.5 * dd$x + rnorm(48)
    r <- run(paste("gaussian y ~ x + (1 | g), no group variance, seed",
                   100 + s),
             frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
    cat("  est sd_g", signif(exp(r$fit$opt$par[length(r$fit$opt$par)]), 3),
        "\n")
  }
}

# ---- 4. two grouping factors that are the same factor -----------------
if (on("dupre")) {
  set.seed(4)
  dd <- data.frame(g = factor(rep(1:15, each = 6)), x = rnorm(90))
  dd$g2 <- factor(as.integer(dd$g) + 100)
  dd$y <- 1 + 0.5 * dd$x + rnorm(15, 0, 1)[dd$g] + rnorm(90)
  run("gaussian y ~ x + (1 | g) + (1 | g2), g2 a relabeled g",
      frm(bf(y ~ x + (1 | g) + (1 | g2)), family = gaussian(), data = dd))
  dd$cnt <- rpois(90, exp(0.5 + 0.3 * dd$x + rnorm(15, 0, 0.7)[dd$g]))
  run("poisson cnt ~ x + (1 | g) + (1 | g2), g2 a relabeled g",
      frm(bf(cnt ~ x + (1 | g) + (1 | g2)), family = poisson(), data = dd))
}

# ---- 5. an RE SD against the residual SD (OLRE) with the check off -----
if (on("olre")) {
  set.seed(5)
  dd <- data.frame(id = factor(1:80), x = rnorm(80))
  dd$y <- 1 + 0.5 * dd$x + rnorm(80)
  run("gaussian y ~ x + (1 | id), one row per id, check_olre = 'ignore'",
      frm(bf(y ~ x + (1 | id)), family = gaussian(), data = dd,
          control = frmtmb_control(check_olre = "ignore")))
  # the same ridge without the structural check: one id per pair of rows
  # with a second, unrelated, singular parameter (a bound)
}

# ---- 6. response on a large scale: rows small in absolute terms --------
# se_empty_row is an absolute 1e-6 in the optimizer's units. A response
# in units of 1e5 gives an identified slope a Hessian row near 1e-8.
if (on("bigscale")) {
  set.seed(6)
  n <- 100
  dd <- data.frame(x = rnorm(n))
  dd$y <- 1e5 * (2 + 0.5 * dd$x + rnorm(n))
  run("y in 1e5 units: y ~ a + b + c * x (a + b a ridge)",
      frm(bf(y ~ a + b + c * x, a ~ 1, b ~ 1, c ~ 1, nl = TRUE),
          family = gaussian(), data = dd,
          start = list(beta = c(1e5, 1e5, 5e4))))
  run("control: y / 1e5, the same model",
      frm(bf(y ~ a + b + c * x, a ~ 1, b ~ 1, c ~ 1, nl = TRUE),
          family = gaussian(), data = transform(dd, y = y / 1e5),
          start = list(beta = c(1, 1, 0.5))))
  run("y in 1e5 units, linear: y ~ x, plus nothing singular (control)",
      frm(bf(y ~ x), family = gaussian(), data = dd))
}
