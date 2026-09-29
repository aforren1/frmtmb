# REVIEW script 09, claim 1 (brms arm): frmtmb.sample's log_lik()
# against brms 2.23.0's log_lik(), cell by cell, at brms's OWN draws.
#
#   Rscript dev/arcovsample-rev-09-brms.R
#
# Own seed (9301) and own designs, chosen to be the ones the worker did
# NOT run: p = 2 and q = 2, weights(), trunc(), cens(), a
# distributional sigma, and set_rescor(TRUE) with a term on each
# response. brms's posterior is irrelevant: every draw is transplanted
# into a frmtmb draws object by parameter name, so both packages report
# a row density at the SAME parameter vector. The transplant is guarded:
# any frmtmb column left unfilled stops the model, and
# posterior_epred() is compared first because it is a function of every
# transplanted parameter.
#
# brms's log_lik undoes its own order(gr, time) sort, so both matrices
# are in the user's row order and are compared without reordering. That
# is asserted rather than assumed, by comparing epred as well.

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
cat("brms ", format(packageVersion("brms")), "  lane build\n", sep = "")

SEED <- 9301L
set.seed(SEED)
lens <- c(5L, 9L, 4L, 8L, 6L, 7L)
dd <- do.call(rbind, lapply(seq_along(lens), function(i) {
  tt <- seq_len(lens[i] + 2L)[-c(2L, 4L)][seq_len(lens[i])]
  data.frame(g = factor(i, levels = seq_along(lens)), t = tt)
}))
n <- nrow(dd)
dd$x <- rnorm(n)
dd$z <- runif(n, -1, 1)
u <- rnorm(length(lens), 0, 0.5)
dd$y <- 0.8 + 0.6 * dd$x + u[as.integer(dd$g)] +
  unlist(lapply(split(seq_len(n), dd$g), function(r) {
    as.numeric(arima.sim(list(ar = c(0.4, 0.2), ma = c(0.3, 0.1)),
                         length(r), sd = 0.5))
  }))
dd$yp <- dd$y + 4
dd$w <- runif(n, 0.4, 2.2)
dd$cc <- sample(c(0L, 0L, 0L, 1L, -1L), n, TRUE)
dd$y2 <- -0.3 + 0.4 * dd$x + rnorm(n, 0, 0.7)
dd <- dd[sample(n), ]
rownames(dd) <- NULL
cat("N = ", n, "  group sizes ", paste(lens, collapse = " "), "  seed ",
    SEED, "\n", sep = "")

# Build the frmtmb draws matrix from a brms draws matrix and a name map.
transplant <- function(fit, bm, map) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  nd <- nrow(bm)
  m <- matrix(NA_real_, nd, length(lab), dimnames = list(NULL, lab))
  for (nm in names(map)) {
    if (!nm %in% lab) stop("map names a column frmtmb does not have: ", nm)
    m[, nm] <- map[[nm]]
  }
  m[, "lp__"] <- 0
  unfilled <- lab[colSums(is.na(m)) > 0L]
  if (length(unfilled)) {
    stop("unfilled frmtmb columns: ", paste(unfilled, collapse = ", "))
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

one <- function(label, bform, bfam, fform, ffam, mapfun, data = dd,
                iter = 300L) {
  cat("\n---- ", label, "\n", sep = "")
  bfit <- suppressWarnings(suppressMessages(
    brms::brm(bform, data = data, family = bfam, chains = 1, iter = iter,
              warmup = 200L, refresh = 0, seed = 4, silent = 2)))
  bm <- as.matrix(brms::as_draws_matrix(bfit))
  fit <- frm(fform, family = ffam, data = data, dry_run = "objective")
  ds <- tryCatch(transplant(fit, bm, mapfun(bm, fit)),
                 error = function(e) { cat("  TRANSPLANT FAILED: ",
                                           conditionMessage(e), "\n",
                                           sep = ""); NULL })
  if (is.null(ds)) return(invisible(NULL))
  ep_b <- brms::posterior_epred(bfit)
  ep_f <- posterior_epred(ds)
  cat("  draws x cols: ", paste(dim(ep_f), collapse = " x "),
      "  epred max|d| ", format(max(abs(ep_f - ep_b)), digits = 6),
      " (epred sd ", format(sd(ep_b), digits = 4), ")\n", sep = "")
  ll_b <- brms::log_lik(bfit)
  ll_f <- log_lik(ds)
  stopifnot(identical(dim(ll_b), dim(ll_f)))
  d <- abs(ll_f - ll_b)
  cat("  log_lik max|d| ", format(max(d), digits = 6),
      "  rel to cell sd ", format(max(d) / sd(ll_b), digits = 4),
      "  bitwise identical: ", identical(ll_f, ll_b), "\n", sep = "")
  # the groups' FIRST rows, where brms's e_s = 0 convention lives
  ord <- order(data$g, data$t)
  first <- ord[!duplicated(data$g[ord])]
  cat("  first row of each group: max|d| ",
      format(max(d[, first]), digits = 6), " over ", length(first),
      " columns\n", sep = "")
  cat("  elpd_loo brms ",
      format(loo::loo(ll_b)$estimates["elpd_loo", "Estimate"],
             digits = 9), "  frmtmb ",
      format(suppressWarnings(loo(ds))$estimates["elpd_loo", "Estimate"],
             digits = 9), "\n", sep = "")
  invisible(NULL)
}

rmap <- function(bm, fit) {
  rn <- grep("^r_g\\[", colnames(bm), value = TRUE)
  out <- list()
  for (nm in rn) out[[nm]] <- bm[, nm]
  out
}

# 1. gaussian arma(2,2): four thetaac, AR before MA
one("gaussian arma(2,2)",
    brms::bf(y ~ x + arma(t, g, p = 2, q = 2)),
    brms::brmsfamily("gaussian"),
    bf(y ~ x + arma(t, g, p = 2, q = 2)), gaussian(),
    function(bm, fit) list(b_Intercept = bm[, "b_Intercept"],
                           b_x = bm[, "b_x"], sigma = bm[, "sigma"],
                           thetaac_1 = bm[, "ar[1]"],
                           thetaac_2 = bm[, "ar[2]"],
                           thetaac_3 = bm[, "ma[1]"],
                           thetaac_4 = bm[, "ma[2]"]))

# 2. student arma(1,1) + (1 | g), nu fitted
one("student arma(1,1) + (1 | g)",
    brms::bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
    brms::brmsfamily("student"),
    bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)), student(),
    function(bm, fit) c(list(b_Intercept = bm[, "b_Intercept"],
                             b_x = bm[, "b_x"], sigma = bm[, "sigma"],
                             nu = bm[, "nu"],
                             theta_1 = log(bm[, "sd_g__Intercept"]),
                             thetaac_1 = bm[, "ar[1]"],
                             thetaac_2 = bm[, "ma[1]"]),
                        rmap(bm, fit)))

# 3. weights() with ma(2)
one("gaussian ma(2), weights(w)",
    brms::bf(y | weights(w) ~ x + ma(t, g, q = 2)),
    brms::brmsfamily("gaussian"),
    bf(y | weights(w) ~ x + ma(t, g, q = 2)), gaussian(),
    function(bm, fit) list(b_Intercept = bm[, "b_Intercept"],
                           b_x = bm[, "b_x"], sigma = bm[, "sigma"],
                           thetaac_1 = bm[, "ma[1]"],
                           thetaac_2 = bm[, "ma[2]"]))

# 4. cens() with arma(2,2)
one("gaussian arma(2,2), cens(cc)",
    brms::bf(y | cens(cc) ~ x + arma(t, g, p = 2, q = 2)),
    brms::brmsfamily("gaussian"),
    bf(y | cens(cc) ~ x + arma(t, g, p = 2, q = 2)), gaussian(),
    function(bm, fit) list(b_Intercept = bm[, "b_Intercept"],
                           b_x = bm[, "b_x"], sigma = bm[, "sigma"],
                           thetaac_1 = bm[, "ar[1]"],
                           thetaac_2 = bm[, "ar[2]"],
                           thetaac_3 = bm[, "ma[1]"],
                           thetaac_4 = bm[, "ma[2]"]))

# 5. trunc() with ar(2)
one("gaussian ar(2), trunc(lb = 0)",
    brms::bf(yp | trunc(lb = 0) ~ x + ar(t, g, p = 2)),
    brms::brmsfamily("gaussian"),
    bf(yp | trunc(lb = 0) ~ x + ar(t, g, p = 2)), gaussian(),
    function(bm, fit) list(b_Intercept = bm[, "b_Intercept"],
                           b_x = bm[, "b_x"], sigma = bm[, "sigma"],
                           thetaac_1 = bm[, "ar[1]"],
                           thetaac_2 = bm[, "ar[2]"]))

# 6. distributional sigma with ar(2)
one("gaussian ar(2), sigma ~ z",
    brms::bf(y ~ x + ar(t, g, p = 2), sigma ~ z),
    brms::brmsfamily("gaussian"),
    bf(y ~ x + ar(t, g, p = 2), sigma ~ z), gaussian(),
    function(bm, fit) list(b_Intercept = bm[, "b_Intercept"],
                           b_x = bm[, "b_x"],
                           b_sigma_Intercept = bm[, "b_sigma_Intercept"],
                           b_sigma_z = bm[, "b_sigma_z"],
                           thetaac_1 = bm[, "ar[1]"],
                           thetaac_2 = bm[, "ar[2]"]))

cat("\n---- labels for the rescor model, to build its map\n")
fr <- frm(bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g)) +
            set_rescor(TRUE), family = gaussian(), data = dd,
          dry_run = "objective")
print(frmtmb::brms_par_labels(fr))

cat("\nDONE\n")
