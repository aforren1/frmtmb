# Reviewer of lane setier: designs that stress se_curvature_real().
#   Rscript dev/setier-rev-probe.R lane|base [part ...]
# Parts:
#   poly   raw polynomials of degree 5..8, x on [1, 2] and on [1e5, 2e5],
#          with and without (1 | g), against lm() and lmer() (seed = degree)
#   curved a curved ridge: y ~ a * b, a ~ 1, b ~ 1 + x + (1 | g); only the
#          products are identified (seeds 1..10)
#   c0k    lane nanse review's c0 + exp(a)^k, a ~ 1 + (1 | g) (seed 77, 3 reps)
#   expb   y ~ a * exp(b * x) with x on [1, 1 + 5e-5]: identified, nearly
#          collinear, where a step of 0.06 SE overflows exp() (seeds 1..5)
#   bound  degree-5 polynomial with an active upper bound on the x5
#          coefficient, with and without (1 | g)
args <- commandArgs(TRUE)
arm <- args[1]
parts <- if (length(args) > 1) args[-1] else
  c("poly", "curved", "c0k", "expb", "bound")
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  merge = c("C:/Users/adf44/source/r/setier-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("arm", arm, find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
catch <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  })
  list(value = v, w = w, m = m)
}
short <- function(x) {
  if (length(x)) substr(paste(x, collapse = " || "), 1, 110) else "none"
}
spec <- function(fit) {
  sdr <- ns$sdr_of(fit)
  an <- fit$cache$se_analysis
  H <- fit$cache$hessian_fixed$H
  if (is.null(H)) H <- tryCatch(fit$obj$he(fit$opt$par), error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(NA)
  D <- sqrt(abs(diag(H))); if (!all(D > 0)) return(NA)
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  min(ev) / max(abs(ev))
}
on <- function(p) p %in% parts

if (on("poly")) {
  cat("\n== poly: raw polynomials\n")
  for (scale in c(1, 1e5)) for (deg in 5:8) for (re in c(FALSE, TRUE)) {
    set.seed(deg)
    d <- data.frame(x = runif(200, 1, 2) * scale, g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x / scale) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    rhs <- paste0("x", 1:deg, collapse = " + ")
    fo <- as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
    r <- tryCatch(catch(frm(fo, data = d)), error = function(e) e)
    if (inherits(r, "error")) {
      cat(sprintf("scale %g deg %d re %s: ERROR %s\n", scale, deg, re,
                  conditionMessage(r)))
      next
    }
    f <- r$value
    fe <- suppressWarnings(fixef(f))
    se <- fe[, "Est.Error"]
    ref <- if (!re) {
      m <- lm(fo, data = d)
      s <- sqrt(diag(vcov(m)) * (200 - length(coef(m))) / 200)
      s[match(rownames(fe), names(s))]
    } else {
      m <- suppressMessages(suppressWarnings(lmer(fo, data = d, REML = FALSE)))
      s <- sqrt(diag(as.matrix(vcov(m))))
      names(s) <- names(fixef(m))
      s[match(rownames(fe), names(s))]
    }
    lost <- ns$sdr_of(f)$se_lost
    rel <- max(abs(se / ref - 1), na.rm = TRUE)
    cat(sprintf(paste0("scale %g deg %d re %-5s code %d | coefs %d (lm/lmer ",
                       "%d) | ev ratio %.2g | lost %s | finite %d | max ",
                       "|SE/ref-1| %.3g | %s\n"),
                scale, deg, re, f$opt$convergence, nrow(fe),
                sum(!is.na(ref)), spec(f),
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "none",
                sum(is.finite(se)), rel, short(c(r$w, r$m))))
  }
}

if (on("curved")) {
  cat("\n== curved ridge y ~ a * b, b ~ 1 + x + (1 | g)\n")
  for (s in 1:10) {
    set.seed(s)
    d <- data.frame(x = rnorm(120), g = factor(rep(1:12, 10)))
    d$y <- 2 * (1 + 0.5 * d$x) + rnorm(12, 0, 0.5)[d$g] + rnorm(120, 0, 0.5)
    r <- catch(frm(bf(y ~ a * b, a ~ 1, b ~ 1 + x + (1 | g), nl = TRUE),
                   data = d))
    f <- r$value
    se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
    nm <- ns$outer_par_names(f)
    ab <- grepl("^(a|b)_", nm)
    lost <- ns$sdr_of(f)$se_lost
    cat(sprintf("seed %2d code %d | ev ratio %.2g | a/b SEs %s | lost %s | %s\n",
                s, f$opt$convergence, spec(f),
                paste(signif(se[ab], 3), collapse = " "),
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "none",
                short(c(r$w, r$m))))
  }
}

if (on("c0k")) {
  cat("\n== c0 + exp(a)^k with (1 | g)\n")
  set.seed(77)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  for (rep in 1:3) {
    r <- catch(frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
                      nl = TRUE), data = dn))
    f <- r$value
    se <- suppressWarnings(fixef(f))[, "Est.Error"]
    lost <- ns$sdr_of(f)$se_lost
    cat(sprintf("rep %d code %d | SE %s | lost %s | %s\n", rep,
                f$opt$convergence, paste(signif(se, 3), collapse = " "),
                paste(names(lost), lost, sep = ":", collapse = ","),
                short(c(r$w, r$m))))
  }
}

if (on("expb")) {
  cat("\n== y ~ a * exp(b * x), x on [1, 1 + 5e-5]\n")
  for (s in 1:5) {
    set.seed(s)
    d <- data.frame(x = 1 + runif(60) * 5e-5)
    d$y <- 2 * exp(0.5 * d$x) + rnorm(60, 0, 1)
    r <- tryCatch(catch(frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                            data = d)),
                  error = function(e) e)
    if (inherits(r, "error")) {
      cat("seed", s, "ERROR", conditionMessage(r), "\n"); next
    }
    f <- r$value
    se <- suppressWarnings(fixef(f))[, "Est.Error"]
    # reference: the exact Hessian inverted as it is
    H <- f$obj$he(f$opt$par)
    ref <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA, 3))
    lost <- ns$sdr_of(f)$se_lost
    cat(sprintf(paste0("seed %d code %d | est %s | ev ratio %.2g | SE %s | ",
                       "solve(H) SE %s | lost %s | %s\n"), s,
                f$opt$convergence, paste(signif(f$opt$par, 4), collapse = " "),
                spec(f), paste(signif(se, 3), collapse = " "),
                paste(signif(ref, 3), collapse = " "),
                paste(names(lost), lost, sep = ":", collapse = ","),
                short(c(r$w, r$m))))
  }
}

if (on("bound")) {
  cat("\n== degree-5 polynomial, upper bound on x5 active\n")
  for (re in c(FALSE, TRUE)) {
    set.seed(5)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
    fo <- as.formula(paste("y ~ x1 + x2 + x3 + x4 + x5", if (re) "+ (1 | g)"))
    f0 <- suppressWarnings(frm(fo, data = d))
    b5 <- fixef(f0)["x5", "Estimate"]
    se5 <- fixef(f0)["x5", "Est.Error"]
    # an upper bound one standard error below the estimate holds it there
    ubv <- b5 - se5
    r <- catch(frm(fo, data = d,
                   prior = set_prior("", class = "b", coef = "x5", ub = ubv)))
    f <- r$value
    fe <- suppressWarnings(fixef(f))
    lost <- ns$sdr_of(f)$se_lost
    # reference: the others with x5 held at the bound
    d$yo <- d$y - ubv * d$x5
    fo2 <- as.formula(paste("yo ~ x1 + x2 + x3 + x4", if (re) "+ (1 | g)"))
    ref <- if (!re) {
      m <- lm(fo2, data = d); sqrt(diag(vcov(m)) * (200 - 5) / 200)
    } else {
      m <- suppressMessages(lmer(fo2, data = d, REML = FALSE))
      sqrt(diag(as.matrix(vcov(m))))
    }
    se <- fe[1:5, "Est.Error"]
    cat(sprintf(paste0("re %s code %d | x5 est %.4g (ub %.4g) | lost %s | SE %s",
                       " | ref (x5 held) %s | %s\n"), re, f$opt$convergence,
                fe["x5", "Estimate"], ubv,
                paste(names(lost), lost, sep = ":", collapse = ","),
                paste(signif(se, 4), collapse = " "),
                paste(signif(ref, 4), collapse = " "), short(c(r$w, r$m))))
  }
}
