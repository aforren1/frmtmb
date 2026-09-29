# Reviewer, claim 1: the Est.Error of a category probability.
#
# The reference is Monte Carlo over the joint covariance the delta method
# uses, with the b positions this script determines FROM THE FRAME rather
# than from smooth_b_idx(), so the reference does not inherit the code's
# own rule:
#   re_formula = NA   -> outer parameters + every smooth/gp/hsgp block's b
#   re_formula = NULL -> outer parameters + every block's b
# At 6000 draws the relative Monte Carlo standard error of one sd is
# 1/sqrt(2 * (N - 1)) = 0.913 percent.
#
#   REVLIB=base Rscript dev/resmooth-rev-estse.R > dev/resmooth-rev-estse-base.txt
#   Rscript dev/resmooth-rev-estse.R > dev/resmooth-rev-estse-lane.txt
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb from:", find.package("frmtmb"), "\n")
NDRAW <- 6000L
cat("MC draws:", NDRAW, "| relative MC se of one sd:",
    sprintf("%.5f", 1 / sqrt(2 * (NDRAW - 1))), "\n\n")

ns <- asNamespace("frmtmb")
gt <- function(nm) get(nm, envir = ns)

# b positions from the frame, independent of the package's own rule
ref_b <- function(fit, keep_re) {
  bl <- fit$frame[["re_blocks"]]
  if (!length(bl)) return(integer(0))
  is_sm <- vapply(bl, function(b) {
    b[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, NA)
  take <- if (keep_re) rep(TRUE, length(bl)) else is_sm
  sort(unique(unlist(lapply(bl[take], `[[`, "b_idx"))))
}

mc_sd <- function(fit, f, idx, ndraw = NDRAW, seed = 11) {
  map <- gt("fit_draw_space")(fit)$map
  v0 <- gt("fit_outer_vector")(fit, map)
  p <- length(v0)
  b0 <- fit$estimates[["b"]]
  V <- if (length(idx)) gt("fd_joint_cov")(fit, map, idx) else
    gt("fit_draw_space")(fit)$V
  if (is.null(V)) return(NULL)
  V <- (V + t(V)) / 2
  L <- t(chol(V + diag(1e-10 * mean(diag(V)) + 0, nrow(V))))
  set.seed(seed)
  m0 <- as.matrix(f(fit))
  out <- matrix(NA_real_, ndraw, length(m0))
  for (i in seq_len(ndraw)) {
    z <- as.vector(L %*% stats::rnorm(nrow(V)))
    g <- gt("fit_set_outer")(fit, v0 + z[seq_len(p)], map)
    if (length(idx)) g$estimates[["b"]][idx] <- b0[idx] + z[p + seq_along(idx)]
    g$cache <- new.env(parent = emptyenv())
    out[i, ] <- as.vector(as.matrix(f(g)))
  }
  list(sd = apply(out, 2, stats::sd), m0 = m0)
}

report <- function(lab, fit, nd, refs = c("NA", "NULL")) {
  bl <- fit$frame[["re_blocks"]]
  cat("==", lab, "\n")
  cat("   blocks:",
      paste(vapply(bl, function(b) b[["covstruct"]], ""), collapse = ","),
      "| n_b:", length(fit$estimates[["b"]]), "\n")
  for (rf in refs) {
    ref <- if (rf == "NA") NA else NULL
    keep <- rf == "NULL"
    ship <- tryCatch(fitted(fit, newdata = nd, re_formula = ref),
                     error = function(e) paste("ERROR:", conditionMessage(e)))
    if (is.character(ship)) { cat("  re_formula =", rf, ship, "\n"); next }
    se_ship <- if (length(dim(ship)) == 3L) {
      as.vector(ship[, "Est.Error", ])
    } else as.vector(ship[, "Est.Error"])
    f <- function(x) {
      gt("fitted_point")(x, nd, ref, "response", NULL, NULL, FALSE)
    }
    idx <- ref_b(fit, keep)
    mc <- tryCatch(mc_sd(fit, f, idx), error = function(e) NULL)
    if (is.null(mc)) { cat("  re_formula =", rf, "MC unavailable\n"); next }
    r <- se_ship / mc$sd
    cat(sprintf("  re_formula = %-4s  b differenced by reference: %d\n",
                rf, length(idx)))
    cat("    shipped Est.Error:", sprintf("%.5f", se_ship), "\n")
    cat("    Monte Carlo sd   :", sprintf("%.5f", mc$sd), "\n")
    cat(sprintf("    ratio ship/MC: min %.4f med %.4f max %.4f\n",
                min(r), stats::median(r), max(r)))
  }
  cat("\n")
}

## ---- A: cumulative, population smooth -------------------------------
set.seed(23)
n <- 240
dA <- data.frame(x = stats::runif(n, -2, 2))
latA <- 1.2 * sin(2 * dA$x) + stats::rlogis(n)
dA$y <- factor(cut(latA, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fA <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                           data = dA))
report("A cumulative y ~ s(x, k = 8)", fA, data.frame(x = c(-1.5, 0, 1.5)))

## ---- B: cumulative, smooth beside (1 | g) ---------------------------
set.seed(31)
ng <- 8; per <- 30
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- 1.0 * sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
report("B cumulative y ~ s(x, k = 8) + (1 | g)", fB,
       dB[c(5, 40, 200), c("x", "g")])

## ---- C: cumulative, t2() --------------------------------------------
set.seed(37)
nC <- 250
dC <- data.frame(x = stats::runif(nC, -2, 2), z = stats::runif(nC, -2, 2))
latC <- sin(dC$x) + 0.6 * dC$z + stats::rlogis(nC)
dC$y <- factor(cut(latC, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fC <- suppressWarnings(frm(bf(y ~ t2(x, z, k = 4)), family = cumulative(),
                           data = dC))
report("C cumulative y ~ t2(x, z, k = 4)", fC,
       data.frame(x = c(-1, 0, 1), z = c(0.5, -0.5, 1)))

## ---- D: cumulative, gp() --------------------------------------------
set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- tryCatch(suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(),
                                    data = dD)),
               error = function(e) {cat("D fit ERROR:", conditionMessage(e),
                                        "\n"); NULL})
if (!is.null(fD)) {
  report("D cumulative y ~ gp(x)", fD, data.frame(x = c(-1.5, 0, 1.5)))
}

## ---- F: cumulative with thres(gr = g) -------------------------------
set.seed(43)
ngF <- 4; perF <- 70
dF <- data.frame(g = factor(rep(seq_len(ngF), each = perF)),
                 x = stats::runif(ngF * perF, -2, 2))
latF <- sin(1.5 * dF$x) + stats::rlogis(nrow(dF))
cutsF <- list(c(-1.5, 0.5), c(-0.5, 1.0), c(-1.0, 0.0), c(-0.2, 1.4))
dF$y <- ordered(vapply(seq_len(nrow(dF)), function(i) {
  cc <- cutsF[[as.integer(dF$g[i])]]
  1L + sum(latF[i] > cc)
}, 1L))
fF <- tryCatch(suppressWarnings(
  frm(bf(y | thres(gr = g) ~ s(x, k = 6)), family = cumulative(), data = dF)),
  error = function(e) {cat("F fit ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(fF)) {
  report("F cumulative y | thres(gr = g) ~ s(x, k = 6)", fF,
         dF[c(3, 80, 150, 250), c("x", "g")])
}

## ---- I: a smooth penalized to near zero -----------------------------
set.seed(47)
nI <- 400
dI <- data.frame(x = stats::runif(nI, -2, 2))
latI <- 0.9 * dI$x + stats::rlogis(nI)      # exactly linear latent mean
dI$y <- factor(cut(latI, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fI <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                           data = dI))
cat("I smooth sd estimates:",
    sprintf("%.6g", exp(fI$estimates[["theta"]])), "\n")
report("I cumulative y ~ s(x, k = 8), near-linear truth", fI,
       data.frame(x = c(-1.5, 0, 1.5)))

## ---- E: distributional sigma ~ s(x), gaussian -----------------------
set.seed(53)
nE <- 300
dE <- data.frame(x = stats::runif(nE, -2, 2))
dE$y <- sin(dE$x) + stats::rnorm(nE, 0, exp(-0.5 + 0.4 * sin(2 * dE$x)))
fE <- suppressWarnings(frm(bf(y ~ s(x, k = 6), sigma ~ s(x, k = 6)),
                           family = gaussian(), data = dE))
ndE <- data.frame(x = c(-1.5, 0, 1.5))
for (rf in c("NA", "NULL")) {
  ref <- if (rf == "NA") NA else NULL
  a <- fitted(fE, newdata = ndE, re_formula = ref)
  cat("E gaussian sigma ~ s(x): fitted() dim",
      paste(dim(a), collapse = "x"), "| re_formula =", rf,
      "| is.matrix:", is.matrix(a), "\n")
  cat("   mu    Est.Error:", sprintf("%.6f", a[, "Est.Error"]), "\n")
  b <- fitted(fE, newdata = ndE, re_formula = ref, dpar = "sigma")
  cat("   sigma Est.Error:", sprintf("%.6f", b[, "Est.Error"]), "\n")
}
# which route? the scalar route is analytic, so no b is differenced
cat("E: matrix route reached?",
    is.matrix(as.matrix(gt("fitted_point")(fE, ndE, NA))), "\n")
for (rf in list(NA, NULL)) {
  f <- function(x) gt("fitted_point")(x, ndE, rf, "response", NULL, NULL,
                                      FALSE)
  mc <- mc_sd(fE, f, ref_b(fE, is.null(rf)))
  a <- fitted(fE, newdata = ndE, re_formula = rf)
  cat("E re_formula =", if (is.null(rf)) "NULL" else "NA",
      "| ratio ship/MC:",
      sprintf("%.4f", a[, "Est.Error"] / mc$sd), "\n")
}
cat("\n")

## ---- H: NO smooth, the identity across builds -----------------------
set.seed(59)
ngH <- 10; perH <- 25
dH <- data.frame(g = factor(rep(seq_len(ngH), each = perH)),
                 x = stats::runif(ngH * perH, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(ngH, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
ndH <- dH[c(2, 60, 180), c("x", "g")]
outH <- list(insample_NA = fitted(fH, re_formula = NA),
             insample_NULL = fitted(fH, re_formula = NULL),
             nd_NA = fitted(fH, newdata = ndH, re_formula = NA),
             nd_NULL = fitted(fH, newdata = ndH, re_formula = NULL))
# a no-smooth gaussian too, for the scalar route
set.seed(61)
dH2 <- data.frame(g = factor(rep(seq_len(ngH), each = perH)),
                  x = stats::runif(ngH * perH, -2, 2))
dH2$y <- 0.8 * dH2$x + stats::rnorm(ngH, 0, 0.7)[dH2$g] +
  stats::rnorm(nrow(dH2), 0, 0.5)
fH2 <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = gaussian(),
                            data = dH2))
outH$g_NA <- fitted(fH2, re_formula = NA)
outH$g_NULL <- fitted(fH2, re_formula = NULL)
saveRDS(outH, file.path("dev", paste0("resmooth-rev-nosmooth-", ARM,
                                      ".rds")))
cat("H no-smooth fitted() saved for", ARM, "\n")
cat("H ordinal Est.Error[1,,1..3]:",
    sprintf("%.10f", outH$insample_NULL[1, "Est.Error", ]), "\n\n")

## ---- G: a nonlinear formula with a smooth in a parameter ------------
set.seed(67)
nG <- 250
dG <- data.frame(x = stats::runif(nG, 0.1, 3), z = stats::runif(nG, -2, 2))
aG <- 2 + 0.5 * sin(2 * dG$z)
dG$y <- aG * exp(-0.8 * dG$x) + stats::rnorm(nG, 0, 0.15)
fG <- tryCatch(suppressWarnings(
  frm(bf(y ~ a * exp(-b * x), a ~ s(z, k = 6), b ~ 1, nl = TRUE),
      family = gaussian(), data = dG)),
  error = function(e) {cat("G fit ERROR:", conditionMessage(e), "\n"); NULL})
if (!is.null(fG)) {
  ndG <- data.frame(x = c(0.5, 1.5), z = c(-1, 1))
  for (rf in list(NA, NULL)) {
    a <- tryCatch(fitted(fG, newdata = ndG, re_formula = rf),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    if (is.character(a)) { cat("G", a, "\n"); next }
    cat("G nl fit, re_formula =", if (is.null(rf)) "NULL" else "NA",
        "| Estimate:", sprintf("%.6f", a[, "Estimate"]),
        "| Est.Error:", sprintf("%.6f", a[, "Est.Error"]), "\n")
    f <- function(x) gt("fitted_point")(x, ndG, rf, "response", NULL, NULL,
                                        FALSE)
    mc <- tryCatch(mc_sd(fG, f, ref_b(fG, is.null(rf)), ndraw = 2000),
                   error = function(e) NULL)
    if (!is.null(mc)) {
      cat("G   MC sd (2000):", sprintf("%.6f", mc$sd), "| ratio:",
          sprintf("%.4f", a[, "Est.Error"] / mc$sd), "\n")
    }
  }
}
cat("\nDONE\n")
