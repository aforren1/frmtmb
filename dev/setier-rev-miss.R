# Reviewer of lane setier: lme4-singular fits of dev/setier-singular.R that
# the lane does not flag (ri20 seed 30, ri6 seed 3, bin15 seed 47, ri20s
# seed 67): the theta row of the finite-difference Hessian, its noise,
# the unit-diagonal spectrum, the at-edge probe, and the SE kept.
#   Rscript dev/setier-rev-miss.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
gen <- function(design, seed) {
  set.seed(seed)
  if (design %in% c("ri6", "ri20", "ri20s")) {
    G <- if (design == "ri6") 6 else 20
    sdg <- if (design == "ri20s") 0.3 else 0
    d <- data.frame(g = factor(rep(seq_len(G), each = 5)), x = rnorm(G * 5))
    d$y <- 1 + 0.5 * d$x + rnorm(G, 0, sdg)[d$g] + rnorm(G * 5)
    list(d = d, ff = y ~ x + (1 | g), fam = gaussian())
  } else {
    d <- data.frame(g = factor(rep(1:15, each = 4)), x = rnorm(60))
    p <- plogis(-0.5 + 0.5 * d$x + rnorm(15, 0, 0)[d$g])
    d$k <- rbinom(60, 10, p)
    d$nk <- 10 - d$k
    list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial())
  }
}
for (cs in list(c("ri20", 30), c("ri6", 3), c("bin15", 47), c("ri20s", 67),
                c("ri20", 1))) {
  g <- gen(cs[1], as.integer(cs[2]))
  m <- character()
  f <- withCallingHandlers(frm(g$ff, family = g$fam, data = g$d),
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  nm <- ns$outer_par_names(f)
  hc <- f$cache$hessian_fixed
  H <- hc$H; E <- hc$E
  j <- match("theta_1", nm)
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat(sprintf("%s seed %s: theta_1 %.3f | row %s | noise row max %.3g | ",
              cs[1], cs[2], f$opt$par[j],
              paste(signif(H[j, ], 3), collapse = " "), max(E[j, ])))
  cat(sprintf("unit-diag ev %s | could_act %s | at_edge(-2) %s | SE theta %.4g | msg %d\n",
              paste(signif(ev, 3), collapse = " "),
              if (exists("se_tier3_could_act", ns))
                ns$se_tier3_could_act(f, H, E, FALSE) else NA,
              if (exists("se_at_edge", ns))
                ns$se_at_edge(f, "theta_1", -1) else NA,
              se[j], length(m)))
  f0 <- f$obj$fn(f$opt$par)
  for (h in c(-2, -5, 2)) {
    q <- f$opt$par; q[j] <- q[j] + h
    cat(sprintf("   theta_1 %+d: nll change %.3g\n", h, f$obj$fn(q) - f0))
  }
}
