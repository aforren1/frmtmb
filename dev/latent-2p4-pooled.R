# Lane `latent`, punch round 3: the coverage POOLED over both blocks.
#
# Round two reported the two blocks side by side and attributed the
# under-coverage to the fresh one. That was the wrong reading: at 200
# replicates the Monte Carlo half-width on a single coefficient is
# about 3 points, so neither block alone can separate a 91 percent row
# from 95. Pooling the 400 can, and it turns a statement about one
# draw into a statement about the DESIGN.
#
# Nothing is re-fitted here. Both files were written by
# dev/latent-2p4-recheck.R at the shipped weight; this reads their raw
# `err_` and `se_` columns, which are already matched to the truth by
# the permutation each row records.
#
#   Rscript dev/latent-2p4-pooled.R

source("dev/latent-env.R")
a <- utils::read.delim("dev/latent-2p4-recheck.tsv", check.names = FALSE)
b <- utils::read.delim("dev/latent-2p4-recheck-oos.tsv",
                       check.names = FALSE)
stopifnot(!any(a$seed %in% b$seed))
cat("tuning block   :", nrow(a), "replicates,", min(a$seed), "to",
    max(a$seed), "\n")
cat("fresh block    :", nrow(b), "replicates,", min(b$seed), "to",
    max(b$seed), "\n")

nm <- c("c2_1", "c2_2", "c2_3", "c3_1", "c3_2", "c3_3",
        "c4_1", "c4_2", "c4_3")
lab <- c("class2:(Intercept)", "class2:x1", "class2:x2",
         "class3:(Intercept)", "class3:x1", "class3:x2",
         "class4:(Intercept)", "class4:x1", "class4:x2")
q <- stats::qnorm(0.975)
row_of <- function(d, k) {
  e <- d[[paste0("err_", nm[k])]]
  s <- d[[paste0("se_", nm[k])]]
  list(cov = abs(e) < q * s, se = s, e = e)
}
ci <- function(x, n) {
  p <- x / n
  h <- q * sqrt(p * (1 - p) / n)
  c(p, p - h, p + h)
}

cat("\n== pooled over all", nrow(a) + nrow(b), "replicates ==\n")
cat(sprintf("  %-20s %6s %9s %8s %18s %8s\n", "coefficient", "n",
            "covered", "rate", "95% interval", "se/sd"))
tot_c <- 0L
tot_n <- 0L
for (k in seq_along(nm)) {
  ra <- row_of(a, k)
  rb <- row_of(b, k)
  cv <- c(ra$cov, rb$cov)
  se <- mean(c(ra$se, rb$se))
  sd <- stats::sd(c(ra$e, rb$e))
  n <- length(cv)
  z <- ci(sum(cv), n)
  tot_c <- tot_c + sum(cv)
  tot_n <- tot_n + n
  flag <- if (z[3L] < 0.95) "  EXCLUDES 95" else ""
  cat(sprintf("  %-20s %6d %9d %7.2f%% %8.2f to %-7.2f %7.4f%s\n",
              lab[k], n, sum(cv), 100 * z[1L], 100 * z[2L], 100 * z[3L],
              se / sd, flag))
}
z <- ci(tot_c, tot_n)
cat(sprintf("\n  %-20s %6d %9d %7.2f%% %8.2f to %-7.2f%s\n", "overall",
            tot_n, tot_c, 100 * z[1L], 100 * z[2L], 100 * z[3L],
            if (z[3L] < 0.95) "  EXCLUDES 95" else ""))

cat("\n== the same rows, per block, for the record ==\n")
cat(sprintf("  %-20s %14s %14s\n", "coefficient", "tuning", "fresh"))
for (k in seq_along(nm)) {
  ra <- row_of(a, k)
  rb <- row_of(b, k)
  cat(sprintf("  %-20s %13.1f%% %13.1f%%\n", lab[k],
              100 * mean(ra$cov), 100 * mean(rb$cov)))
}

cat("\n== is the shortfall a tail, or the whole distribution? ==\n")
cat("  z = err / se. sd(z) is 1 when the reported standard error\n")
cat("  matches the spread; IQR(z)/1.349 is the same thing robustly.\n")
cat(sprintf("  %-20s %8s %8s %8s %10s\n", "coefficient", "sd(z)",
            "robust", "ratio", "kurtosis"))
for (k in c(3L, 6L, 9L)) {
  ra <- row_of(a, k)
  rb <- row_of(b, k)
  zz <- c(ra$e, rb$e) / c(ra$se, rb$se)
  kur <- mean((zz - mean(zz))^4) / stats::sd(zz)^4
  rob <- stats::IQR(zz) / 1.349
  cat(sprintf("  %-20s %8.4f %8.4f %8.4f %10.3f\n", lab[k],
              stats::sd(zz), rob, stats::sd(zz) / rob, kur))
}
