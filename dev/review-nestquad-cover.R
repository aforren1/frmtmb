# Fix round 1, M1: coverage of the Wald and bootstrap intervals that
# remain on nested quadrature fits. Bernoulli, 15 schools x 3 classes x 4
# rows, beta = (0.2, 0.7), sds 0.8 and 0.8.
#   Rscript dev/review-nestquad-cover.R wald <reps>
#   Rscript dev/review-nestquad-cover.R boot <reps> <nsim>
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
a <- commandArgs(TRUE)
what <- a[1]; reps <- as.integer(a[2]); nsim <- as.integer(a[3] %||% "0")
truth <- c(0.2, 0.7, log(0.8), log(0.8))
sim <- function(seed) {
  set.seed(seed)
  d <- expand.grid(r = 1:4, cls = 1:3, school = 1:15)
  d$school <- factor(d$school)
  d$cls_uni <- factor(paste0(d$school, "_", d$cls))
  d$x <- rnorm(nrow(d))
  eta <- truth[1] + truth[2] * d$x + rnorm(15, 0, 0.8)[d$school] +
    rnorm(45, 0, 0.8)[d$cls_uni]
  d$y <- rbinom(nrow(d), 1, plogis(eta))
  d
}
hit <- matrix(NA, reps, 4); nwarn <- 0; nerr <- 0
t0 <- proc.time()[[3]]
for (i in seq_len(reps)) {
  d <- sim(1000 + i)
  f <- try_msg(frm(bf(y ~ x + (1 | school/cls_uni)) + bernoulli(), data = d,
                   quadrature = !grepl("lap", what)))
  if (inherits(f, "error")) { nerr <- nerr + 1; next }
  if (length(attr(f, "warns"))) nwarn <- nwarn + 1
  ci <- if (what %in% c("wald", "waldlap")) {
    try(confint(f), silent = TRUE)
  } else {
    try(suppressWarnings(confint(f, method = "boot", nsim = nsim, seed = i)),
        silent = TRUE)
  }
  if (inherits(ci, "try-error")) { nerr <- nerr + 1; next }
  ci <- ci[seq_len(4), 1:2, drop = FALSE]
  hit[i, ] <- ci[, 1] <= truth & truth <= ci[, 2]
}
cat(sprintf("COVER %s reps %d nsim %d errors %d fits-with-warnings %d time %.0fs\n",
            what, reps, nsim, nerr, nwarn, proc.time()[[3]] - t0))
cat("COVER rows:", paste(c("Intercept", "x", "theta_school", "theta_class"),
                         sprintf("%.3f (n %d)", colMeans(hit, na.rm = TRUE),
                                 colSums(!is.na(hit))), collapse = "; "), "\n")
