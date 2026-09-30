# Reviewer, lane ordinal: frmtmb.sample draws on the new structures. For
# each shape, posterior_epred() at a draw (which goes through the inverse
# map to the internal vector) against brms's R-side density at that
# draw's STORED columns (b_Intercept[...], b_x, disc), so a wrong inverse
# shows. Data seed 20261009, sampler seed 5.
# Output: dev/ordinal-rev-log-draws.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(20261009)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$y2 <- 1L + (stats::rlogis(n) + d$x > 0) + (stats::rlogis(n) + d$x > 1)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
shapes <- list(
  cum_stz = list(f = y ~ x, fam = cumulative(threshold = "sum_to_zero"),
                 dens = "cumulative"),
  acat_stz_probit = list(f = y ~ x, fam = acat("probit", threshold = "sum_to_zero"),
                         dens = "acat", link = "probit"),
  cratio_equi = list(f = y ~ x, fam = cratio(threshold = "equidistant"),
                     dens = "cratio"),
  sratio_equi_gr = list(f = y | thres(gr = h) ~ x,
                        fam = sratio(threshold = "equidistant"), dens = "sratio",
                        gr = TRUE),
  cum_stz_gr = list(f = y | thres(gr = h) ~ x,
                    fam = cumulative(threshold = "sum_to_zero"), dens = "cumulative",
                    gr = TRUE),
  cum_equi_disc = list(f = bf(y ~ x, disc ~ 0 + z),
                       fam = cumulative(threshold = "equidistant"),
                       dens = "cumulative", disc = TRUE)
)
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  cat("\n==", nm, "\n")
  tryCatch({
    fit <- frm(s$f, family = s$fam, data = d)
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
    v <- variables(ds)
    cat("  variables:", v[!grepl("^(lp__|lprior)", v)], "\n")
    M <- as.matrix(ds)
    ep <- posterior_epred(ds)
    worst <- 0
    for (i in c(1L, 57L, nrow(M))) {
      eta <- d$x * M[i, "b_x"]
      disc <- if (isTRUE(s$disc)) exp(d$z * M[i, "b_disc_z"]) else 1
      lk <- s$link %||% "logit"
      if (isTRUE(s$gr)) {
        P <- matrix(NA_real_, n, 5)
        for (g in c("p", "q")) {
          cols <- grep(paste0("^b_Intercept\\[", g, ","), colnames(M))
          rows <- which(d$h == g)
          th <- M[i, cols]
          P[rows, seq_len(length(th) + 1L)] <- get(paste0("d", s$dens),
            asNamespace("brms"))(seq_len(length(th) + 1L), eta = eta[rows],
            thres = matrix(th, length(rows), length(th), byrow = TRUE),
            disc = disc, link = lk)
        }
      } else {
        th <- M[i, grep("^b_Intercept\\[", colnames(M))]
        P <- get(paste0("d", s$dens), asNamespace("brms"))(1:5, eta = eta,
          thres = matrix(th, n, length(th), byrow = TRUE), disc = disc, link = lk)
      }
      worst <- max(worst, max(abs(ep[i, , ] - P), na.rm = TRUE))
      if (i == 1L && grepl("stz", nm)) {
        cat("  draw 1 threshold sums:", if (isTRUE(s$gr)) {
          c(sum(M[1, grep("^b_Intercept\\[p,", colnames(M))]),
            sum(M[1, grep("^b_Intercept\\[q,", colnames(M))]))
        } else sum(th), "\n")
      }
      if (i == 1L && grepl("equi", nm)) {
        dl <- grep("^delta", colnames(M), value = TRUE)
        cat("  draw 1 delta columns:", dl, "=", M[1, dl], "\n")
      }
    }
    cat("  max |posterior_epred - brms at stored columns| over 3 draws:", worst, "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
cat("\n== mv: equidistant y, sum_to_zero y2 ==\n")
tryCatch({
  fit <- frm(bf(y ~ x) + bf(y2 ~ x), data = d,
             family = list(cumulative(threshold = "equidistant"),
                           sratio(threshold = "sum_to_zero")))
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
  v <- variables(ds)
  cat("  variables:", v[!grepl("^(lp__|lprior)", v)], "\n")
  M <- as.matrix(ds)
  th <- M[1, grep("^b_y_Intercept\\[", colnames(M))]
  cat("  y equidistant check, diff(th) - delta_y:", diff(th) - M[1, "delta_y"], "\n")
  th2 <- M[1, grep("^b_y2_Intercept\\[", colnames(M))]
  cat("  y2 sum:", sum(th2), "\n")
}, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
