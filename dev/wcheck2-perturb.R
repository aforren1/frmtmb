# Lane wcheck2: the new test block's fits under rounding. The data of
# wu_data() in tests/testthat/test-weights-unbounded.R (seed 20261020),
# with x and the response scaled by 1 + 2e-15 k, k = 0..24. For each
# family: the guard-absent fit (the point mass at weight 3, the other
# rows' excess -2) must end at code 0 with no warning, and the refusal
# (the mass at 0.5, the other rows' excess +4) must hold on every draw.
#   Rscript dev/wcheck2-perturb.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
wu_data <- function(seed = 20261020) {
  set.seed(seed)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  tt = factor(rep(seq_len(nt), ng)),
                  x = rnorm(n), xs = rnorm(n), w = runif(n, 0.5, 2))
  d$obs <- factor(seq_len(n))
  eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$ly <- log(d$y)
  d$k <- rpois(n, exp(eta))
  d$lsig <- rep(log(0.7), n)
  d$fac <- factor(rep(c("a", "b"), length.out = n))
  d
}
d0 <- wu_data()
n <- nrow(d0)
i0 <- seq_len(n) %% 5 == 0
i1 <- seq_len(n) %% 7 == 0 & !i0
cases <- list(
  hurdle_lognormal = list(fam = hurdle_lognormal(), y = "yz", mass = i0),
  hurdle_gamma = list(fam = hurdle_gamma(), y = "yz", mass = i0),
  zero_inflated_beta = list(fam = zero_inflated_beta(), y = "yb",
                            mass = i0),
  zero_one_inflated_beta = list(fam = zero_one_inflated_beta(), y = "yo",
                                mass = i0 | i1))
res <- list()
for (k in 0:24) {
  s <- 1 + 2e-15 * k
  d <- d0
  d$x <- d0$x * s
  d$yz <- ifelse(i0, 0, d0$y * s)
  yb <- stats::plogis(d0$ly - mean(d0$ly)) * s
  d$yb <- ifelse(i0, 0, yb)
  d$yo <- d$yb
  d$yo[i1] <- 1
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    m <- cs$mass
    fo <- as.formula(paste(cs$y, "| weights(w) ~ x + (1 | obs)"))
    d$w <- ifelse(m, 3, 1 - 2 / sum(!m))
    o <- outcome(frm(fo, family = cs$fam, data = d))
    tag <- if (startsWith(o, "fit")) {
      sprintf("code %s, %s warnings",
              sub(".* conv ([0-9]+) .*", "\\1", o),
              sub(".* warn ([0-9]+).*", "\\1", o))
    } else substr(o, 1, 40)
    res[[paste(nm, "absent")]] <- c(res[[paste(nm, "absent")]], tag)
    d$w <- ifelse(m, 0.5, 1 + 4 / sum(!m))
    o <- outcome(frm(fo, family = cs$fam, data = d))
    res[[paste(nm, "refused")]] <- c(res[[paste(nm, "refused")]],
                                     if (startsWith(o, "REFUSED")) "REFUSED"
                                     else substr(o, 1, 40))
  }
}
for (nm in names(res)) {
  tb <- table(res[[nm]])
  cat(sprintf("%-36s %s\n", nm,
              paste(sprintf("%s x%d", names(tb), as.integer(tb)),
                    collapse = "; ")))
}
cat("DONE\n")
