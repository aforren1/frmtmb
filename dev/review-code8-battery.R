# The code-8 battery of the bigx review (dev/reviews/2026-11-06-bigx.md,
# confirmation): products and powers of shifted columns, raw polynomials,
# and y ~ w far from zero, each judged against the same model on centered
# or orthogonal columns.
#   CODE8_ARM=before|after Rscript dev/code8-battery.R <label> [parts]
# parts: comma list of prod, poly, yw, ywx (yw extended: seeds 911-940
# and 1e10). Output dev/code8-out/battery-<label>.tsv
a <- commandArgs(TRUE)
lab <- a[1]
parts <- if (length(a) > 1) strsplit(a[2], ",")[[1]] else
  c("prod", "poly", "yw")
source("dev/code8-common.R")
rows <- list()
rec <- function(case, r, ref_ll) {
  si <- stop_info(r)
  cv <- r$wc[grepl("frmtmb_conv_code", r$wc)]
  rows[[length(rows) + 1L]] <<- cbind(
    data.frame(case = case), si,
    dll = if (r$err) NA else ll_of(r$f) - ref_ll,
    w_code = sum(grepl("frmtmb_conv_code", r$wc)),
    m_atmax = sum(grepl("frmtmb_converged_at_max", r$mc)),
    w_other = sum(!grepl("frmtmb_conv_code", r$wc)),
    conds = cond_txt(r))
}
if ("prod" %in% parts) {
  for (s in 1:5) for (m in c(1e4, 1e5, 1e6)) {
    set.seed(40 + s); n <- 300
    z1 <- rnorm(n); z2 <- rnorm(n)
    d <- data.frame(z1 = z1, z2 = z2, w1 = z1 + m, w2 = z2 + m)
    d$y <- 0.5 + 0.4 * z1 - 0.3 * z2 + 0.3 * z1 * z2 + rnorm(n, 0, 0.5)
    d$yp <- rpois(n, exp(0.2 + 0.3 * z1 + 0.05 * z1^2))
    fz <- quiet(frm(bf(y ~ z1 * z2), data = d))
    rec(sprintf("gauss w1*w2 m=%g s%d", m, 40 + s),
        rfit(bf(y ~ w1 * w2), data = d), ll_of(fz))
    fz <- quiet(frm(bf(yp ~ z1 + I(z1^2)), family = poisson(), data = d))
    rec(sprintf("pois w1+I(w1^2) m=%g s%d", m, 40 + s),
        rfit(bf(yp ~ w1 + I(w1^2)), family = poisson(), data = d),
        ll_of(fz))
  }
}
if ("poly" %in% parts) {
  for (s in 501:503) {
    set.seed(s); n <- 200
    xx <- runif(n, 1000, 1010)
    d <- data.frame(xx = xx, yy = sin((xx - 1000) / 2) + rnorm(n, 0, 0.3))
    for (deg in c(5, 6, 8)) {
      fr <- as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))
      fo <- as.formula(sprintf("yy ~ poly(xx, %d)", deg))
      fz <- quiet(frm(bf(fo), data = d))
      rec(sprintf("poly deg %d s%d", deg, s), rfit(bf(fr), data = d),
          ll_of(fz))
    }
  }
}
yw <- function(seeds, ms) {
  for (s in seeds) for (m in ms) for (fam in c("gaussian", "poisson")) {
    set.seed(s); n <- 150; z <- rnorm(n)
    y <- if (fam == "gaussian") 0.5 + 0.4 * z + rnorm(n, 0, 0.5) else
      rpois(n, exp(0.5 + 0.4 * z))
    d <- data.frame(y = y, z = z, w = z + m)
    F <- if (fam == "gaussian") gaussian() else poisson()
    fz <- quiet(frm(bf(y ~ z), family = F, data = d))
    rec(sprintf("%s y~w m=%g s%d", fam, m, s),
        rfit(bf(y ~ w), family = F, data = d), ll_of(fz))
  }
}
if ("yw" %in% parts) yw(901:910, c(1e7, 1e8, 1e9))
if ("ywx" %in% parts) {
  yw(911:940, c(1e7, 1e8, 1e9))
  yw(901:920, 1e10)
}
X <- do.call(rbind, rows)
dir.create("dev/review-code8-out", showWarnings = FALSE)
write.table(X, sprintf("dev/review-code8-out/battery-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
c8 <- X[!is.na(X$code) & X$code != 0, ]
cat("fits", nrow(X), "nonzero code", nrow(c8), "of which false conv",
    sum(grepl("false", c8$msg)), "\n")
print(table(code = X$code, useNA = "ifany"))
if (nrow(c8)) {
  cat("gain range", format(range(c8$gain, na.rm = TRUE), digits = 3),
      " NA gains", sum(is.na(c8$gain)),
      " dll range", format(range(c8$dll, na.rm = TRUE), digits = 3),
      " se not finite", sum(!c8$se_fin), "\n")
  print(c8[, c("case", "code", "gain", "se_fin", "n_lost", "dll", "w_code",
               "m_atmax", "w_other")], row.names = FALSE)
}
cat("DONE\n")
