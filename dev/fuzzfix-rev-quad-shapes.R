# Reviewer: which multi-block shapes does quadrature = TRUE accept,
# refuse, or leave half Laplace? One arm per process (some are costly).
#   Rscript dev/fuzzfix-rev-quad-shapes.R <lib|base> <arm>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
arm <- a[2]

set.seed(31)
ng <- if (length(a) > 2) as.integer(a[3]) else 8L; m <- 3
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
n <- nrow(d)
d$x <- rnorm(n)
d$cell <- interaction(d$ga, d$gb, drop = TRUE)
d$cellu <- factor(paste0("c", as.integer(d$cell)))   # unique labels
d$tt <- factor(rep(1:3, length.out = n))             # crossed with ga
d$f4 <- factor(as.integer(d$cell) %% 4)              # constant in a cell
eta <- 0.3 + 0.5 * d$x + rnorm(ng, 0, 0.8)[d$ga] +
  rnorm(nlevels(d$cell), 0, 0.6)[d$cell]
d$y <- rbinom(n, 1, plogis(eta))
d$k <- rpois(n, exp(eta))
d$yg <- eta + rnorm(n, 0, 0.7)
# one gb level per ga: ga:gb has as many levels as ga
d$gb1 <- factor(rep("i", n))
# partially crossed: one inner level spans two outer levels
d$gbx <- as.character(d$cellu); d$gbx[d$ga == "2" & d$gb == "i"] <- "c1"
d$gbx <- factor(d$gbx)

B <- bernoulli()
arms <- list(
  nested = list(bf(y ~ x + (1 | ga/gb)), B),
  explicit = list(bf(y ~ x + (1 | ga) + (1 | ga:gb)), B),
  unique_labels = list(bf(y ~ x + (1 | ga) + (1 | cellu)), B),
  equal_levels = list(bf(y ~ x + (1 | ga/gb1)), B),
  inner_first_written = list(bf(y ~ x + (1 | cellu) + (1 | ga)), B),
  by_split = list(bf(y ~ x + (1 | gr(ga, by = gb))), B),
  by_split_inner = list(bf(y ~ x + (1 | ga) + (1 | gr(cellu, by = gb))), B),
  by_split_inner4 = list(bf(y ~ x + (1 | ga) + (1 | gr(cellu, by = f4))), B),
  mu_sigma = list(bf(yg ~ x + (1 | ga), sigma ~ (1 | ga)), gaussian()),
  mv = list(mvbf(bf(y ~ x + (1 | ga)) + B, bf(k ~ x + (1 | ga)) + poisson()),
            NULL),
  mv_nested = list(mvbf(bf(y ~ x + (1 | ga/gb)) + B,
                        bf(k ~ x + (1 | ga/gb)) + poisson()),
                   NULL),
  crossed = list(bf(y ~ x + (1 | ga) + (1 | tt)), B),
  partial_cross = list(bf(y ~ x + (1 | ga) + (1 | gbx)), B),
  mu_sigma_nested = list(bf(yg ~ x + (1 | ga/gb), sigma ~ (1 | ga)),
                         gaussian())
)
ar <- arms[[arm]]
fam <- ar[[2]]
t0 <- proc.time()[[3]]
ws <- character()
fq <- tryCatch(withCallingHandlers(
  if (is.list(fam) && !inherits(fam, "family")) {
    frm(ar[[1]], family = fam, data = d, quadrature = TRUE)
  } else frm(ar[[1]], family = fam, data = d, quadrature = TRUE),
  warning = function(w) { ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning") },
  message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
el <- proc.time()[[3]] - t0
if (inherits(fq, "error")) {
  cat(sprintf("ARM %s REFUSED/ERROR (%.1fs): %s\n", arm, el,
              conditionMessage(fq)))
} else {
  fl <- suppressWarnings(suppressMessages(
    frm(ar[[1]], family = fam, data = d)))
  sm <- tryCatch(suppressWarnings(suppressMessages(summary(fq))),
                 error = function(e) e)
  nb <- vapply(fq$frame$re_blocks, function(b) b$n_levels, 1)
  cat(sprintf(paste0("ARM %s ok (%.1fs) blocks %s left_random %d ",
                     "logLik quad %.8f laplace %.8f diff %.3e summary %s ",
                     "warnings %d\n"), arm, el,
              paste(nb, collapse = "+"), length(fq$obj$env$random),
              as.numeric(logLik(fq)), as.numeric(logLik(fl)),
              as.numeric(logLik(fq)) - as.numeric(logLik(fl)),
              if (inherits(sm, "error")) conditionMessage(sm) else "ok",
              length(ws)))
  if (length(ws)) cat("  warnings:", unique(substr(ws, 1, 100)), sep = "\n  ")
}
