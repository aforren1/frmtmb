# Lane wcheck2, backlog 4998: asym_laplace "72 of 80 right censored"
# with (1 | g4). The data of dev/wcheck-cases.R part cens (seed
# 20261021, n = 80). That script's loop leaves d$ca and d$yac at its LAST
# value, k = 80, so the (1 | g4) arm fits 80 of 80 rows right censored.
#   Rscript dev/wcheck2-4998.R <base|lane|lib>
source("dev/wcheck2-common.R")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$ya <- 1 + 0.5 * d$x + rnorm(n)
for (k in c(72, 80)) {
  o <- order(d$ya, decreasing = TRUE)
  ca <- rep("none", n)
  ca[o[seq_len(k)]] <- "right"
  qa <- d$ya[o[k]]
  d$ca <- ca
  d$yac <- pmin(d$ya, qa)
  cat(sprintf("k = %d: %d right censored at %.4f\n", k, sum(ca == "right"),
              qa))
  for (fo in c("yac | cens(ca) ~ x + (1 | g4)", "yac | cens(ca) ~ x",
               "yac | cens(ca) ~ x + (1 | obs)")) {
    cat(sprintf("  asym_laplace %-32s %s\n", fo,
                outcome(frm(as.formula(fo), family = asym_laplace(),
                            data = d))))
    cat(sprintf("  gaussian     %-32s %s\n", fo,
                outcome(frm(as.formula(fo), family = gaussian(), data = d))))
  }
}
cat("DONE\n")
