# Reviewer, claim 3 recovery: the reviewer's rerun of every replicate
# (dev/fams2-rev-coverage.R into dev/fams2-rev-cov/) against the
# worker's dev/fams2-cov/, file by file, and the coverage recomputed
# from the rerun with Wilson intervals.
w <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-cov"
r <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-cov"
fw <- sort(list.files(w, pattern = "[.]rds$"))
fr <- sort(list.files(r, pattern = "[.]rds$"))
cat("files: worker", length(fw), " reviewer", length(fr),
    " same names:", identical(fw, fr), "\n")
same <- 0; maxd <- 0; nerr <- 0; nwarn <- 0
for (f in intersect(fw, fr)) {
  a <- readRDS(file.path(w, f)); b <- readRDS(file.path(r, f))
  if (!is.null(b$error)) nerr <- nerr + 1
  if (length(b$warn)) nwarn <- nwarn + 1
  if (identical(a, b)) same <- same + 1 else if (is.null(a$error) && is.null(b$error)) {
    maxd <- max(maxd, abs(a$est - b$est), abs(a$lwr - b$lwr), abs(a$upr - b$upr))
  }
}
cat("identical replicate files:", same, " max |difference| where not:", maxd,
    "\nreviewer errors:", nerr, " with a warning:", nwarn, "\n")
wilson <- function(k, n, z = stats::qnorm(0.975)) {
  p <- k / n
  c <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(c - h, c + h)
}
res <- lapply(file.path(r, fr), readRDS)
fam <- vapply(res, `[[`, "", "family")
for (f in unique(fam)) {
  rr <- res[fam == f]
  ok <- vapply(rr, function(x) is.null(x$error), NA)
  cov <- do.call(rbind, lapply(rr[ok], `[[`, "cover"))
  est <- do.call(rbind, lapply(rr[ok], `[[`, "est"))
  truth <- rr[[which(ok)[1]]]$truth
  cat(sprintf("\n%s: %d fits\n", f, sum(ok)))
  for (p in colnames(cov)) {
    k <- sum(cov[, p]); n <- nrow(cov); wi <- wilson(k, n)
    bias <- mean(est[, p]) - truth[[p]]
    cat(sprintf("  %-18s cover %.3f [%.3f, %.3f]%s  bias/mcse %+.2f  bias/sd %+.3f\n",
                p, k / n, wi[1], wi[2], if (wi[1] > 0.95 || wi[2] < 0.95) " MISS" else "",
                bias / (sd(est[, p]) / sqrt(n)), bias / sd(est[, p])))
  }
}
