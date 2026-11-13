# Lane gpstop: the cost of the relinearized start, in ONE process on the
# lane library, the lane against the mutant "norelin"
# (gp_body_relin_applies() returns FALSE: the 0.81.0 route), interleaved
# row by row with the order alternating, so load hits both arms. Rows: the
# two 96-row sets of dev/gpstop-gpfresh.R. Args: set (fresh | lane96).
a <- commandArgs(TRUE)
set <- a[1]
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_body_relin_applies
set_arm <- function(arm) {
  rv_put("gp_body_relin_applies",
         if (arm == "norelin") function(...) FALSE else orig)
}
if (set == "fresh") {
  source(file.path(GS_WT, "dev/review-fitedge-designs.R"))
  rows <- expand.grid(seed = 1:2, n = c(60L, 120L), ci = 1:3,
                      design = rv_designs, stringsAsFactors = FALSE)
  getd <- function(r) rv_mk(r$design, r$seed, r$n, r$ci)
  getsp <- function(r) rv_spell[[r$design]]
  lab <- function(r) c("sin", "wig", "exp")[r$ci]
} else {
  curves <- c("sin", "wig", "exp", "step")
  mk <- function(seed, n, ci, body) {
    set.seed(18000 + 100 * seed + n + ci)
    x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
    f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
                exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
    y <- switch(body,
      A = f + 0.3 * z + 0.1 * rnorm(n),
      B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
      C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n),
      D = 1 + f * w + 0.1 * rnorm(n))
    data.frame(x, z, w, y)
  }
  spell <- list(
    A = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             fam = gaussian()),
    B = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             fam = poisson()),
    C = list(nl = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             fam = gaussian()),
    D = list(nl = bf(y ~ b + a * w, a ~ gp(x), b ~ 1, nl = TRUE),
             fam = gaussian()))
  rows <- expand.grid(seed = 1:3, n = c(60L, 120L), ci = 1:4,
                      design = names(spell), stringsAsFactors = FALSE)
  getd <- function(r) mk(r$seed, r$n, r$ci, r$design)
  getsp <- function(r) spell[[r$design]]
  lab <- function(r) curves[r$ci]
}
tot <- c(lane = 0, norelin = 0)
for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  d <- getd(r)
  sp <- getsp(r)
  arms <- if (i %% 2) c("lane", "norelin") else c("norelin", "lane")
  res <- list()
  for (arm in arms) {
    set_arm(arm)
    t0 <- proc.time()[[3]]
    q <- fe_fit(sp$nl, data = d, family = sp$fam)
    el <- proc.time()[[3]] - t0
    tot[arm] <- tot[arm] + el
    ll <- if (inherits(q$fit, "error")) NA else as.numeric(logLik(q$fit))
    res[[arm]] <- c(ll = ll, sec = el,
                    runs = if (inherits(q$fit, "error") ||
                                 is.null(q$fit$opt$gp_second)) 0 else
                      q$fit$opt$gp_second[["runs"]])
  }
  cat(sprintf(paste0("BCOST %s %s %s %d %d lane ll %.6f sec %.2f runs %d |",
                     " norelin ll %.6f sec %.2f runs %d | diff %.6f\n"),
              set, r$design, lab(r), r$n, r$seed, res$lane[["ll"]],
              res$lane[["sec"]], res$lane[["runs"]], res$norelin[["ll"]],
              res$norelin[["sec"]], res$norelin[["runs"]],
              res$lane[["ll"]] - res$norelin[["ll"]]))
  flush.console()
}
set_arm("lane")
cat(sprintf("BCOSTSUM %s rows %d lane %.1f s norelin %.1f s ratio %.2f\n",
            set, nrow(rows), tot[["lane"]], tot[["norelin"]],
            tot[["lane"]] / tot[["norelin"]]))
cat("DONE\n")
