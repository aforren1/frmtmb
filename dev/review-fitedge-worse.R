# Review of lane fitedge, item 1: the rows where the lane reads below
# rellib-r18, reproduced. Args: arm (base | lane). On the lane, each row is
# also fitted with gp_body_start() off (0.79.1's start on the lane's build), so
# "keep the better of two descents" can be read with its cost, in one process.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
source(file.path(FE_WT, "dev/review-fitedge-designs.R"))
ns <- asNamespace("frmtmb")
gbs <- if (arm == "lane") get("gp_body_start", envir = ns)
put <- function(v) {
  if (arm != "lane") return(invisible(NULL))
  unlockBinding("gp_body_start", ns)
  assign("gp_body_start", v, envir = ns)
  lockBinding("gp_body_start", ns)
}
curves <- c("sin", "wig", "exp", "step")
mk96 <- function(seed, n, ci, body) {
  set.seed(18000 + 100 * seed + n + ci)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
  f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  y <- switch(body,
    B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
    C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n))
  data.frame(x, z, w, y)
}
fB <- bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)
fC <- bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE)
rows <- list(
  list("lane96 B step 60 1", fB, poisson(), function() mk96(1, 60L, 4L, "B")),
  list("lane96 B wig 120 3", fB, poisson(), function() mk96(3, 120L, 2L, "B")),
  list("lane96 C exp 60 2", fC, gaussian(), function() mk96(2, 60L, 3L, "C")),
  list("fresh K exp 60 1", rv_spell$K$nl, rv_spell$K$fam,
       function() rv_mk("K", 1, 60L, 3L)),
  list("fresh G exp 60 1", rv_spell$G$nl, rv_spell$G$fam,
       function() rv_mk("G", 1, 60L, 3L)),
  list("fresh J wig 120 2", rv_spell$J$nl, rv_spell$J$fam,
       function() rv_mk("J", 2, 120L, 2L)),
  list("fresh J exp 120 2", rv_spell$J$nl, rv_spell$J$fam,
       function() rv_mk("J", 2, 120L, 3L)))
for (r in rows) {
  d <- r[[4]]()
  arms <- if (arm == "lane") c("on", "off") else "base"
  res <- list()
  for (a1 in arms) {
    if (a1 == "off") put(function(...) NULL)
    t0 <- proc.time()[[3]]
    q <- fe_fit(r[[2]], data = d, family = r[[3]])
    tt <- proc.time()[[3]] - t0
    put(gbs)
    ll <- as.numeric(logLik(q$fit))
    res[[a1]] <- c(ll = ll, t = tt)
    cat(sprintf("WR %s | %s %s ll %.6f conv %d nwarn %d %.2f s | theta %s\n",
                r[[1]], arm, a1, ll, q$fit$opt$convergence, length(q$warn), tt,
                paste(round(q$fit$opt$par[names(q$fit$opt$par) == "theta"], 3),
                      collapse = " ")))
  }
  if (arm == "lane") {
    cat(sprintf("WR %s | best of two %.6f, cost %.2f s against %.2f s (%.2f x)\n",
                r[[1]], max(res$on["ll"], res$off["ll"]),
                res$on["t"] + res$off["t"], res$on["t"],
                (res$on["t"] + res$off["t"]) / res$on["t"]))
  }
}
cat("DONE\n")
