# Lane wcheck2, backlog 4968: hurdle and zero-inflated continuous
# families with weights and (1 | obs). For each family and weight arm:
# the check's verdict (frm()) and the Laplace objective's profile along
# the spike dpar's log (every other parameter re-optimized), with the
# slope per unit toward the limit (-log sigma, or log phi / log shape).
# The zero (and one) rows are a point mass whose own effect neither
# grows nor curves, so the predicted slope is rate * E_nz, where E_nz is
# the excess of the other rows alone (rate 1 for sigma, 1/2 for the
# beta's phi and the gamma's shape).
#   Rscript dev/wcheck2-4968.R <base|lane|lib> <fam> [quad]
# fam: hl (hurdle_lognormal), hg (hurdle_gamma), zib (zero_inflated_beta),
# zoib (zero_one_inflated_beta). Data of dev/reviews/wcheck-review-misc.R
# part quad (seed 20261111, n = 80).
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
fam_key <- a[2]
quad <- identical(a[3], "quad")
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g = factor(rep(1:16, each = 5)), t = rep(1:5, 16))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$yz <- ifelse(runif(n) < 0.2, 0, exp(d$y / 2))
d$yzb <- plogis(d$y - 1)
d$yzb[1:10] <- 0
d$yzo <- d$yzb
d$yzo[11:16] <- 1
spec <- switch(fam_key,
  hl = list(name = "hurdle_lognormal", fam = hurdle_lognormal(), y = "yz",
            zero = c(0), toward_zero = TRUE, rate = 1),
  hg = list(name = "hurdle_gamma", fam = hurdle_gamma(), y = "yz",
            zero = c(0), toward_zero = FALSE, rate = 0.5),
  zib = list(name = "zero_inflated_beta", fam = zero_inflated_beta(),
             y = "yzb", zero = c(0), toward_zero = FALSE, rate = 0.5),
  zoib = list(name = "zero_one_inflated_beta",
              fam = zero_one_inflated_beta(), y = "yzo", zero = c(0, 1),
              toward_zero = FALSE, rate = 0.5))
yv <- d[[spec$y]]
pm <- yv %in% spec$zero
n_nz <- sum(!pm)
cat(sprintf("%s: %d point-mass rows, %d other rows%s\n", spec$name,
            sum(pm), n_nz, if (quad) ", quadrature" else ""))
arms <- list(
  list(lab = "all w 1.5", w = rep(1.5, n)),
  list(lab = "other rows E +4, point mass w 0.5",
       w = ifelse(pm, 0.5, 1 + 4 / n_nz)),
  list(lab = "other rows E -2, point mass w 3",
       w = ifelse(pm, 3, 1 - 2 / n_nz)),
  list(lab = "other rows E -2, point mass w 1",
       w = ifelse(pm, 1, 1 - 2 / n_nz)))
g <- if (spec$toward_zero) -c(2, 4, 6, 8, 10, 12) else c(2, 4, 6, 8, 10, 12)
for (arm in arms) {
  d$w <- arm$w
  E_nz <- sum(arm$w[!pm] - 1)
  E_all <- sum(arm$w - 1)
  fo <- as.formula(paste(spec$y, "| weights(w) ~ x + (1 | obs)"))
  cat(sprintf("--- %s: E_nz %.4g, sum(w) - n %.4g, predicted slope %.4g\n",
              arm$lab, E_nz, E_all, spec$rate * E_nz))
  cat("  frm():", outcome(frm(fo, family = spec$fam, data = d,
                               quadrature = quad)), "\n")
  if (quad) next
  ob <- q(frm(fo, family = spec$fam, data = d, dry_run = "objective"))
  v <- profile_disp(ob, g)
  # slope per unit toward the limit
  gg <- if (spec$toward_zero) -g else g
  show_profile(sprintf("  profile along %s", if (spec$toward_zero)
    "-log sigma" else "log of the spike dpar"), gg, v)
}
cat("DONE\n")
