# Lane illse2, item 4: an active upper bound on the top coefficient of a
# raw polynomial (as dev/setier-rev-bound.R, at degrees 6 to 8). The
# constrained maximum is lm() without the top column and with it as an
# offset at the bound, so the other coefficients' SEs (ML scale) are
# lm()'s from that fit, through its QR.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-p4bound.R <degs> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
degs <- eval(parse(text = a[1]))
seeds <- eval(parse(text = a[2]))
for (deg in degs) for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
  xn <- paste0("x", 1:deg)
  for (k in 1:deg) d[[xn[k]]] <- d$x^k
  top <- xn[deg]
  fo <- as.formula(paste("y ~", paste0(xn, collapse = " + ")))
  l0 <- lm(fo, data = d)
  b <- coef(l0)[[top]]
  s <- sqrt(diag(vcov(l0)))[[top]]
  ub <- b - s
  r <- cap(frm(fo, data = d,
               prior = set_prior("", class = "b", coef = top, ub = ub)))
  if (is_err(r$v)) { cat(deg, seed, "ERROR", r$v, "\n"); next }
  f <- r$v
  d$off <- ub * d[[top]]
  fo1 <- as.formula(paste("y ~", paste0(xn[-deg], collapse = " + "),
                          "+ offset(off)"))
  l1 <- lm(fo1, data = d)
  n <- nrow(d)
  k <- length(coef(l1))
  se_ref <- sqrt(diag(vcov(l1)) * (n - k) / n)
  fe <- suppressWarnings(fixef(f))
  se <- fe[sub("^[(]Intercept[)]$", "Intercept", names(se_ref)),
           "Est.Error"]
  lost <- ns$sdr_of(f)$se_lost
  cat(sprintf(paste0("deg %d seed %d | code %d | %s at bound %s | ",
                     "logLik - lm %+.2e | lost %s | root %s | ",
                     "SE / lm %s | warn %s\n"),
              deg, seed, f$opt$convergence, top,
              isTRUE(abs(fe[top, "Estimate"] - ub) <= 1e-8 * abs(ub)),
              as.numeric(logLik(f)) - as.numeric(logLik(l1)),
              if (length(lost)) paste(names(lost), lost, collapse = ",")
              else "none",
              !is.null(ns$sdr_of(f)$cov_fixed_root),
              ratio_str(unname(se / se_ref)), wtag(r$w)))
}
