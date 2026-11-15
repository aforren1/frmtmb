# Short check of lane silent2's fix round, F3 and 5214 (1).
#   Rscript dev/reviews/silent2-short-f3.R
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "\n")
ws <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(val = val, w = w)
}
show <- function(r) {
  for (cnd in r$w) {
    cat("WARNING class:", paste(class(cnd), collapse = ","), "\n  ",
        gsub("\n", " ", conditionMessage(cnd)), "\n")
  }
}
nan_pars <- function(fit) {
  fx <- suppressWarnings(fixef(fit))
  rownames(fx)[is.nan(fx[, "Est.Error"])]
}

cat("## F3: cs_sr_acat seed 1\n")
gen <- function(seed, n = 300, K = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}
r <- ws(frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()),
            data = gen(1)))
cat("logLik", format(as.numeric(logLik(r$val)), digits = 10), "\n")
show(r)
cat("fixef NaN Est.Error:", nan_pars(r$val), "\n")

cat("## three-component test of test-ordinal-mixture.R\n")
set.seed(20261005 + 11)
n <- 400
x <- rnorm(n)
rnorm(n)
sample(c("a", "b"), n, TRUE)
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
d3 <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))), x = x)
r3 <- ws(frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                         cumulative()),
             data = d3, control = frmtmb_control(starts = 0)))
cat("logLik", format(as.numeric(logLik(r3$val)), digits = 10), "\n")
show(r3)
cat("fixef NaN Est.Error:", nan_pars(r3$val), "\n")
fx <- suppressWarnings(fixef(r3$val))
print(signif(fx[, 1:2], 4))

cat("## 5214 (1)\n")
omx_flat_data <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
dre <- transform(omx_flat_data(9), g = factor(rep_len(1:20, 300)))
fam <- mixture(cumulative(), sratio())
says <- function(r) {
  any(grepl("frmtmb_control(starts = 10)",
            vapply(r$w, conditionMessage, ""), fixed = TRUE))
}
degen <- function(r) {
  any(grepl("degenerate boundary", vapply(r$w, conditionMessage, ""),
            fixed = TRUE))
}
a <- ws(frm(bf(y ~ x + z + (1 | g)), family = fam, data = dre))
cat("RE, starts unset:  degenerate", degen(a), " names starts", says(a), "\n")
b <- ws(frm(bf(y ~ x + z + (1 | g)), family = fam, data = dre,
            control = frmtmb_control(starts = 0)))
cat("RE, starts = 0:    degenerate", degen(b), " names starts", says(b), "\n")
c0 <- ws(frm(bf(y ~ x + z), family = fam, data = omx_flat_data(16)))
cat("no RE, default:    degenerate", degen(c0), " names starts", says(c0), "\n")
c1 <- ws(frm(bf(y ~ x + z), family = fam, data = omx_flat_data(16),
             control = frmtmb_control(starts = 0)))
cat("no RE, starts = 0: degenerate", degen(c1), " names starts", says(c1), "\n")
cat("DONE\n")
