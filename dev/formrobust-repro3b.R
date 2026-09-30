# Item 3 on draws: posterior_predict(newdata = ) with NA responses under
# ar(cov = FALSE), frmtmb.sample. Seed 31 data, sampler seed 3.
LIB <- Sys.getenv("FORMROBUST_LIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/stan-cache")
suppressMessages(library(frmtmb.sample))
cat("frmtmb.sample from", find.package("frmtmb.sample"), "\n")
set.seed(31)
G <- 30; Tn <- 8
d <- expand.grid(t = 1:Tn, g = factor(1:G))
d$x <- rnorm(nrow(d))
e <- as.vector(apply(matrix(rnorm(G * Tn), Tn, G), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$y <- 1 + 0.5 * d$x + e
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ x + ar(t, g, p = 1)), family = gaussian(), data = d,
             chains = 1, iter = 600, refresh = 0, seed = 3)))
nd <- d[d$g %in% c("1", "2"), ]
nd$y[nd$g == "1" & nd$t >= 5] <- NA
r <- tryCatch({
  set.seed(2)
  posterior_predict(ds, newdata = nd)
}, error = function(e) conditionMessage(e))
if (is.character(r)) {
  cat("posterior_predict: ERROR:", r, "\n")
} else {
  i5 <- which(nd$g == "1" & nd$t == 5)
  i6 <- which(nd$g == "1" & nd$t == 6)
  i4 <- which(nd$g == "1" & nd$t == 4)
  cat("dim", dim(r), " anyNA", anyNA(r), "\n")
  sdr <- apply(r, 2, sd)
  cat("sd row 6 / sd row 5:", format(sdr[i6] / sdr[i5], digits = 6),
      " sd row 4 / sd row 5:", format(sdr[i4] / sdr[i5], digits = 6), "\n")
  ep <- posterior_epred(ds, newdata = nd)
  cat("epred anyNA:", anyNA(ep), "\n")
}
