# Reviewer helper: add the fix-round-1 cases to dev/cdfs-rev-brmsfit-r1.R
f <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/cdfs-rev-brmsfit-r1.R"
s <- readLines(f)
i <- grep("^cases\\$Gamma_trunc <- list\\(", s)
add <- c(
  "# fix round 1: exgaussian, right above q85, at sigma / beta about 2 and",
  "# about 25 (the regime RTMBdist::dexgauss() approximated)",
  "for (ub in c(2, 25)) {",
  "  y <- 1 + 0.4 * x + rnorm(n, 0, 0.5) + rexp(n, ub / 0.5)",
  "  d <- data.frame(x = x, y = y, cc = 0, y2 = y)",
  "  q85 <- quantile(y, 0.85)",
  "  d$cc[y > q85] <- 1; d$y[d$cc == 1] <- q85",
  "  cases[[paste0(\"exgaussian_u\", ub)]] <- list(d = d, f = frmtmb::exgaussian(),",
  "                                              b = brms::exgaussian())",
  "}",
  "# fix round 1: poisson, interval bins inside trunc(lb = 1, ub = 20)",
  "y <- rpois(3 * n, exp(1.3 + 0.3 * rep(x, 3)))",
  "xx <- rep(x, 3)",
  "keep <- which(y >= 1 & y <= 20)[1:n]",
  "d <- data.frame(x = xx[keep], y = y[keep], cc = 0)",
  "d$y2 <- d$y",
  "ib <- sample(n, 100)",
  "d$cc[ib] <- 2",
  "d$y[ib] <- pmax(d$y2[ib] - 2, 1)",
  "d$cc[d$cc == 2 & d$y2 <= d$y] <- 0",
  "cases$poisson_trunc_interval <- list(d = d, f = poisson(), b = poisson(),",
  "                                     discrete = TRUE, tr2 = TRUE)")
s <- append(s, add, after = i - 1)
j <- grep("^  lhs <- if \\(isTRUE\\(cs\\$trunc\\)\\)", s)
s[j] <- paste0("  lhs <- if (isTRUE(cs$tr2)) \"y | cens(cc, y2) + trunc(lb = 1, ub = 20)\" else\n",
               "    if (isTRUE(cs$trunc)) \"y | trunc(lb = 0.5, ub = 6)\" else")
writeLines(s, f)
