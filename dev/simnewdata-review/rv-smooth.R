# Is the smooth fix complete? For each smooth construction: which blocks
# simulate(re_formula = NA) redraws, which blocks predict(re_formula = NA)
# keeps, and the empirical check that the draws centre on
# fitted(re_formula = NA) with the spread of what is redrawn.
#   Rscript dev/simnewdata-review/rv-smooth.R > .../log/smooth[-base].txt
source("dev/simnewdata-review/rv-prelude.R")
set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)
d$y2 <- d$y + 0.5 * rnorm(n)
d$ys <- sin(2 * pi * d$x) + rnorm(n, 0, exp(-1 + 1.5 * d$x))

nsim <- 1500
probe <- function(label, f, dpar = "mu") {
  fit <- tryCatch(suppressWarnings(frm(f, data = d)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-26s FIT ERROR: %s\n", label,
                substr(conditionMessage(fit), 1, 150)))
    return(invisible())
  }
  bl <- fit$frame$re_blocks
  cs <- vapply(bl, `[[`, "", "covstruct")
  kept <- integer(0)
  for (lp in fit$frame$linpreds) {
    ed <- frmtmb:::lp_eta_design(fit, lp, NULL, FALSE, FALSE)
    kept <- c(kept, match(vapply(ed$sm_blocks, function(b) b$b_idx[1], 0),
                          vapply(bl, function(b) b$b_idx[1], 0)))
  }
  redraw <- if (!rv_base) frmtmb:::sim_re_plan(fit, NA)$blocks else NA
  s <- tryCatch(as.matrix(simulate(fit, nsim = nsim, seed = 2,
                                   re_formula = NA)),
                error = function(e) e)
  if (inherits(s, "error")) {
    cat(sprintf("%-26s SIM ERROR: %s\n", label,
                substr(conditionMessage(s), 1, 150)))
    return(invisible())
  }
  mu <- fitted(fit, re_formula = NA)[, "Estimate"]
  sdr <- apply(s, 1, sd)
  z <- (rowMeans(s) - mu) / (sdr / sqrt(nsim))
  sg <- as.vector(frm_linpred(fit, dpar = "sigma", type = "response",
                              re_formula = NA))
  if (length(sg) == 1L) sg <- rep(sg, nrow(s))
  cat(sprintf("%-26s blocks %s | pred(NA) keeps %s | sim(NA) redraws %s | mean max|z| %.1f rms|z| %.2f | rowsd/sigma median %.3f\n",
              label, paste(seq_along(cs), cs, sep = ":", collapse = ","),
              paste(sort(unique(kept)), collapse = ","),
              paste(redraw, collapse = ","), max(abs(z)),
              sqrt(mean(z^2)), stats::median(sdr / sg)))
}
probe("s(x)", bf(y ~ s(x)))
probe("s(x, by = f)", bf(y ~ f + s(x, by = f)))
probe("t2(x, z)", bf(y ~ t2(x, z)))
probe("te(x, z)", bf(y ~ te(x, z)))
probe("s(g, bs = re)", bf(y ~ s(x) + s(g, bs = "re")))
probe("s(x, g, bs = fs)", bf(y ~ s(x, g, bs = "fs", k = 5)))
probe("gp(x)", bf(y ~ gp(x)))
probe("gp(x, by = f)", bf(y ~ f + gp(x, by = f)))
probe("gp(x, k = 10) hsgp", bf(y ~ gp(x, k = 10, c = 5 / 4)))
probe("sigma ~ s(x)", bf(ys ~ s(x), sigma ~ s(x)))
probe("s(x) + (1 | g)", bf(y ~ s(x) + (1 | g)))
probe("t2(x, g, re)", bf(y ~ t2(x, g, bs = c("tp", "re"))))
cat("multivariate:",
    try_msg(simulate(suppressWarnings(frm(bf(mvbind(y, y2) ~ s(x)),
                                          data = d)),
                     nsim = 2, re_formula = NA)), "\n")
cat("multivariate bf+bf:",
    try_msg(simulate(suppressWarnings(frm(bf(y ~ s(x)) + bf(y2 ~ s(x)),
                                          data = d)),
                     nsim = 2, re_formula = NA)), "\n")
