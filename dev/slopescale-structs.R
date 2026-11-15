# Lane slopescale, item 1018, brief step 2: a random slope on a column
# of small spread, with no fixed effect on it, across covariance
# structures. Each case is fitted on x * s (s = 1e-4) and on x (spread
# 1); the two are the same model in other units, so logLik must agree,
# the slope's VarCorr sd and its ranef()/coef() must be the spread-1
# values over s, the fitted values and predict() draws (one seed) must
# agree.
# Usage: Rscript slopescale-structs.R before|after [seeds]
a <- commandArgs(TRUE)
arm <- a[1]
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:3
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "arm", arm, "\n")
S <- 1e-4

mk <- function(seed) {
  set.seed(seed)
  ng <- 20; np <- 12
  g <- factor(rep(seq_len(ng), each = np))
  h <- factor(sample(seq_len(ng), ng * np, TRUE))
  a <- factor(rep(1:5, each = ng * np / 5))
  tf <- factor(rep(1:4, length.out = ng * np))
  x0 <- rnorm(ng * np)
  u0 <- rnorm(ng, 0, 0.5)[g]
  u1 <- rnorm(ng, 0, 0.6)[g]
  ua <- rnorm(5, 0, 0.4)[a]
  ut <- rnorm(ng * 4, 0, 0.5)[(as.integer(g) - 1) * 4 + as.integer(tf)]
  eta <- 0.5 + u0 + u1 * x0 + 0.5 * ua * x0 + ut * x0
  sig <- exp(-0.3 + rnorm(ng, 0, 0.4)[g] * x0)
  y <- rnorm(length(eta), eta, sig)
  data.frame(y = y, g = g, h = h, a = a, tf = tf, x = x0)
}
A <- local({
  set.seed(99)
  M <- matrix(rnorm(20 * 30), 20)
  K <- crossprod(t(M)) / 30
  dimnames(K) <- list(as.character(1:20), as.character(1:20))
  K
})

cases <- list(
  us      = list(f = bf(y ~ 1 + (1 + x | g))),
  us0     = list(f = bf(y ~ 1 + (0 + x | g))),
  dbar    = list(f = bf(y ~ 1 + (1 + x || g))),
  cs      = list(f = bf(y ~ 1 + cs(1 + x | g))),
  ar1     = list(f = bf(y ~ 1 + (1 | g) + ar1(0 + tf:x | g))),
  homdiag = list(f = bf(y ~ 1 + (1 | g) + homdiag(0 + x | g))),
  grnocor = list(f = bf(y ~ 1 + (1 + x | gr(g, cor = FALSE)))),
  grcov   = list(f = bf(y ~ 1 + (1 + x | gr(g, cov = A))),
                 data2 = list(A = A)),
  student = list(f = bf(y ~ 1 + (1 + x | gr(g, dist = "student")))),
  mm      = list(f = bf(y ~ 1 + (1 + x | mm(g, h)))),
  nested  = list(f = bf(y ~ 1 + (1 + x | a/g))),
  dpar    = list(f = bf(y ~ 1 + (1 | g), sigma ~ 1 + (0 + x | g))),
  idlink  = list(f = bf(y ~ 1 + (1 + x | p | g), sigma ~ 1 + (1 | p | g)))
)

cnt <- new.env()
run <- function(f, d, data2) {
  cnt$w <- character(0)
  cnt$m <- character(0)
  r <- try(withCallingHandlers(
    frm(f, family = gaussian(), data = d, data2 = data2),
    warning = function(w) {
      cnt$w <- c(cnt$w, class(w)[1])
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cnt$m <- c(cnt$m, class(m)[1])
      invokeRestart("muffleMessage")
    }), silent = TRUE)
  list(fit = r, w = cnt$w, m = cnt$m)
}
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
# sds named by group and coefficient; a slope's sd is in units of 1/s
sds <- function(f, col = "Estimate") {
  v <- VarCorr(f)
  unlist(lapply(names(v), function(gn) {
    sdm <- v[[gn]]$sd
    if (is.null(sdm)) return(NULL)
    stats::setNames(sdm[, col], paste0(gn, ":", rownames(sdm)))
  }))
}
# standard errors: the fixed effects' (none is on the slope's column)
# and the sds' Est.Error, the slope's in units of 1/s
ses <- function(f1, fs) {
  c1 <- summary(f1)$coefficients
  cs_ <- summary(fs)$coefficients
  fe <- max(unlist(lapply(names(c1), function(dp) {
    abs(cs_[[dp]][, "Std. Error"] / c1[[dp]][, "Std. Error"] - 1)
  })))
  e1 <- sds(f1, "Est.Error")
  es <- sds(fs, "Est.Error")
  es[slope(names(es))] <- es[slope(names(es))] * S
  ok <- is.finite(e1) & is.finite(es)
  c(fe = fe, vc = if (any(ok)) max(abs(es[ok] / e1[ok] - 1)) else NA,
    nan_mismatch = sum(is.finite(e1) != is.finite(es)))
}
slope <- function(nm) grepl("x$", nm)

rows <- list()
for (cn in names(cases)) for (seed in seeds) {
  cs <- cases[[cn]]
  d1 <- mk(seed)
  ds <- d1
  ds$x <- d1$x * S
  r1 <- run(cs$f, d1, cs$data2)
  rs <- run(cs$f, ds, cs$data2)
  if (inherits(r1$fit, "try-error") || inherits(rs$fit, "try-error")) {
    cat(cn, seed, "ERROR",
        if (inherits(r1$fit, "try-error")) r1$fit else rs$fit, "\n")
    rows[[length(rows) + 1L]] <- data.frame(case = cn, seed = seed,
      short = NA, conv = NA, nwarn = NA, boundary = NA, vc = NA,
      re = NA, fit = NA, pred = NA, se_fe = NA, se_vc = NA, se_nan = NA)
    next
  }
  f1 <- r1$fit
  fs <- rs$fit
  short <- as.numeric(logLik(f1)) - as.numeric(logLik(fs))
  s1 <- sds(f1)
  ss <- sds(fs)
  ss[slope(names(ss))] <- ss[slope(names(ss))] * S
  vc <- tryCatch(rel(ss, s1[names(ss)]), error = function(e) NA)
  # ranef() and coef(): one matrix (or data frame) per grouping factor,
  # a column per coefficient; the slope columns are in units of 1/s
  recmp <- function(l1, ls) {
    max(vapply(names(l1), function(gn) {
      a1 <- as.matrix(l1[[gn]])
      as_ <- as.matrix(ls[[gn]])
      sc <- ifelse(slope(colnames(as_)), S, 1)
      as_ <- sweep(as_, 2, sc, "*")
      max(abs(as_ - a1)) / max(abs(a1))
    }, 0))
  }
  re <- tryCatch(max(recmp(ranef(f1), ranef(fs)), recmp(coef(f1), coef(fs))),
                 error = function(e) NA)
  ft <- tryCatch(max(abs(fitted(fs)[, "Estimate"] -
                           fitted(f1)[, "Estimate"])) / stats::sd(d1$y),
                 error = function(e) NA)
  pr <- tryCatch({
    set.seed(7)
    p1 <- predict(f1)[, "Estimate"]
    set.seed(7)
    ps <- predict(fs)[, "Estimate"]
    max(abs(ps - p1)) / stats::sd(d1$y)
  }, error = function(e) NA)
  se <- tryCatch(ses(f1, fs), error = function(e) c(fe = NA, vc = NA,
                                                     nan_mismatch = NA))
  rows[[length(rows) + 1L]] <- data.frame(
    case = cn, seed = seed, short = short, conv = fs$opt$convergence,
    nwarn = length(rs$w),
    boundary = "frmtmb_boundary_fit" %in% rs$m, vc = vc, re = re,
    fit = ft, pred = pr, se_fe = se[["fe"]], se_vc = se[["vc"]],
    se_nan = se[["nan_mismatch"]])
}
r <- do.call(rbind, rows)
options(width = 160)
print(r, digits = 4)
cat("\nshort > 1e-6:", sum(r$short > 1e-6, na.rm = TRUE), "of", nrow(r),
    "; with no warning:", sum(r$short > 1e-6 & r$nwarn == 0, na.rm = TRUE),
    "; errors:", sum(is.na(r$short)), "; max short",
    max(r$short, na.rm = TRUE), "\n")
nn <- r[r$case != "nested", ]
cat("non-nested max: vc", signif(max(nn$vc, na.rm = TRUE), 3),
    "re", signif(max(nn$re, na.rm = TRUE), 3),
    "fit", signif(max(nn$fit, na.rm = TRUE), 3),
    "pred", signif(max(nn$pred, na.rm = TRUE), 3),
    "se_fe", signif(max(nn$se_fe, na.rm = TRUE), 3),
    "se_vc", signif(max(nn$se_vc, na.rm = TRUE), 3),
    "se_nan mismatches", sum(nn$se_nan, na.rm = TRUE), "\n")
