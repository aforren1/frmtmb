# REVIEW script 10: edge constructions that could break the new path,
# plus the NaN seen in script 04's synthetic laplace case.
#
#   Rscript dev/arcovsample-rev-10-edges.R

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
`%||%` <- function(a, b) if (is.null(a)) b else a

re_prior_of <- function(sh) {
  tot <- 0
  for (bk in sh$frame[["re_blocks"]] %||% list()) {
    f <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll
    tot <- tot + as.numeric(f(sh$estimates[["b"]][bk$b_idx],
                              sh$estimates[["theta"]][bk$theta_idx], bk))
  }
  tot
}

mkd <- function(fit, vals, nd = 3L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(NA_real_, nd, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) {
    if (!nm %in% lab) stop("no such column: ", nm)
    m[, nm] <- vals[[nm]]
  }
  m[, "lp__"] <- 0
  un <- lab[colSums(is.na(m)) > 0L]
  if (length(un)) stop("unfilled: ", paste(un, collapse = ", "))
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

check <- function(lab, ds) {
  ll <- log_lik(ds)
  fit <- frmtmb.sample:::draws_base_fit(ds)
  idx <- frmtmb.sample:::draws_par_index(fit)
  r <- numeric(nrow(ds$draws))
  for (k in seq_along(r)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    nll <- as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
    r[k] <- sum(ll[k, ]) - (-nll - re_prior_of(sh))
  }
  cat("  ", lab, ": dim ", paste(dim(ll), collapse = "x"),
      "  max|rowsum - objective| ", format(max(abs(r)), digits = 8),
      "  rel ", format(max(abs(r) / abs(rowSums(ll))), digits = 4),
      "  finite ", all(is.finite(ll)), "\n", sep = "")
}

cat("--- A. a group with ONE row, and a group with two\n")
set.seed(501L)
dA <- rbind(
  data.frame(g = "a", t = 1L),
  data.frame(g = "b", t = 1:2),
  data.frame(g = "c", t = c(1L, 3L, 4L, 7L)),
  data.frame(g = "d", t = 1:9))
dA$g <- factor(dA$g)
dA$x <- rnorm(nrow(dA))
dA$y <- 0.4 + 0.5 * dA$x + rnorm(nrow(dA), 0, 0.7)
dA <- dA[sample(nrow(dA)), ]
for (pq in list(c(1L, 0L), c(2L, 0L), c(0L, 2L), c(2L, 2L))) {
  p <- pq[1L]; q <- pq[2L]
  form <- if (q == 0L) bf(y ~ x + ar(t, g, p = p)) else
    if (p == 0L) bf(y ~ x + ma(t, g, q = q)) else
      bf(y ~ x + arma(t, g, p = p, q = q))
  f <- frm(form, family = gaussian(), data = dA, dry_run = "objective")
  v <- list(b_Intercept = 0.4, b_x = 0.5, sigma = 0.7)
  for (j in seq_len(p + q)) v[[paste0("thetaac_", j)]] <- 0.3 - 0.05 * j
  check(paste0("p=", p, " q=", q), mkd(f, v))
}

cat("\n--- B. two responses, ar() on ONE of them, NO rescor\n")
dB <- dA
dB$y2 <- -0.2 + 0.3 * dB$x + rnorm(nrow(dB), 0, 0.6)
f <- frm(bf(y ~ x + ar(t, g, p = 2)) + bf(y2 ~ x), family = gaussian(),
         data = dB, dry_run = "objective")
cat("  arma_cond_resp: ", paste(frmtmb::arma_cond_resp(f),
                                collapse = ","), "\n", sep = "")
print(frmtmb::brms_par_labels(f))
check("mv, no rescor", mkd(f, list(b_y_Intercept = 0.4, b_y_x = 0.5,
                                   sigma_y = 0.7, b_y2_Intercept = -0.2,
                                   b_y2_x = 0.3, sigma_y2 = 0.6,
                                   thetaac_1 = 0.35, thetaac_2 = 0.15)))

cat("\n--- C. ar() on BOTH responses, no rescor\n")
f2 <- frm(bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g)),
          family = gaussian(), data = dB, dry_run = "objective")
cat("  arma_cond_resp: ", paste(frmtmb::arma_cond_resp(f2),
                                collapse = ","), "\n", sep = "")
check("mv, both cond", mkd(f2, list(b_y_Intercept = 0.4, b_y_x = 0.5,
                                    sigma_y = 0.7, b_y2_Intercept = -0.2,
                                    b_y2_x = 0.3, sigma_y2 = 0.6,
                                    thetaac_1 = 0.35, thetaac_2 = 0.2)))

cat("\n--- D. one response cov = FALSE, the other cov = TRUE\n")
f3 <- frm(bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g, cov = TRUE)),
          family = gaussian(), data = dB, dry_run = "objective")
cat("  arma_cond_resp: ", paste(frmtmb::arma_cond_resp(f3),
                                collapse = ","), "\n", sep = "")
lab <- c(frmtmb::brms_par_labels(f3), "lp__")
fd <- structure(list(stanfit = NULL,
                     draws = matrix(0, 2L, length(lab),
                                    dimnames = list(NULL, lab)),
                     fit = f3), class = "frmtmb_draws")
cat("  log_lik: ", conditionMessage(tryCatch(log_lik(fd),
                                             error = identity)), "\n",
    sep = "")

cat("\n--- E. a REAL laplace draws object, both builds' concern\n")
dE <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dE$x <- rnorm(nrow(dE))
dE$y <- 0.5 + 0.4 * dE$x + rnorm(nrow(dE), 0, 0.8)
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
dsl <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x + ar(t, g) + (1 | g)), family = gaussian(),
             data = dE, chains = 1, iter = 300, refresh = 0, seed = 3,
             laplace = TRUE)))
cat("  is laplace: ", frmtmb.sample:::draws_is_laplace(dsl), "\n",
    sep = "")
set.seed(9)
pp <- tryCatch(posterior_predict(dsl), error = identity)
if (inherits(pp, "condition")) {
  cat("  posterior_predict: ", class(pp)[1L], ": ",
      conditionMessage(pp), "\n", sep = "")
} else {
  cat("  posterior_predict dim ", paste(dim(pp), collapse = "x"),
      "  NaN cells ", sum(!is.finite(pp)), " of ", length(pp), "\n",
      sep = "")
}
cat("  log_lik: ",
    conditionMessage(tryCatch({ log_lik(dsl); simpleError("NO ERROR") },
                              error = identity)), "\n", sep = "")
cat("\nDONE\n")
