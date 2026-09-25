# Reviewer 2, item 6: summarise the base-build coverage run from the
# RDS files alone. Usage: Rscript r2-cov-summary.R <lib>
a <- commandArgs(trailingOnly = TRUE)
lib <- if (length(a)) a[[1]] else "base"
dir <- file.path("dev/phase3b-review2", paste0("cov-", lib))
# optional seed range, <lo> <hi>, and the file-time cutoff: a file
# written at or after 18:30 on 2026-09-24 is void (battery loss, 18:39)
fs <- list.files(dir, "^seed-[0-9]+[.]rds$", full.names = TRUE)
sn <- as.integer(gsub("[^0-9]", "", basename(fs)))
if (length(a) >= 3) fs <- fs[sn >= as.integer(a[[2]]) & sn <= as.integer(a[[3]])]
mt <- file.mtime(fs)
cut <- as.POSIXct("2026-09-24 18:30:00")
cat("files written at or after 18:30 (void):", sum(mt >= cut), "; latest",
    format(max(mt)), "\n")
fs <- fs[mt < cut]
R <- lapply(fs, function(f) tryCatch(readRDS(f), error = function(e)
  list(seed = NA, error = paste("unreadable:", f))))
cat("unreadable:", sum(vapply(R, function(r) is.na(r$seed), TRUE)), "\n")
seeds <- vapply(R, function(r) r$seed, 0)
err <- vapply(R, function(r) !is.null(r$error), TRUE)
cat(sprintf("%s: %d files, seeds %s; errors %d\n", lib, length(fs),
            paste(range(seeds), collapse = " to "), sum(err)))
gaps <- setdiff(seq(min(seeds), max(seeds)), seeds)
cat("gaps in the seed range:", if (length(gaps)) paste(gaps, collapse = " ") else "none", "\n")
R <- R[!err]
cat("build:", unique(vapply(R, function(r) r$eam_path, "")), unique(vapply(R, function(r) r$eam_ver, "")), "\n")
cat("convergence code 0:", sum(vapply(R, function(r) r$conv == 0, TRUE)), "of", length(R), "\n")
wilson <- function(k, n) {
  z <- 1.959964; p <- k / n
  c(p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)),
    p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n)
}
get <- function(row, col) vapply(R, function(r) r$ci[row, col], 0)
n <- length(R)
row <- function(lab, nm, truth) {
  est <- get(nm, "est"); lwr <- get(nm, "lwr"); upr <- get(nm, "upr")
  se <- (upr - lwr) / (2 * qnorm(0.975))
  k <- sum(lwr < truth & truth < upr); w <- wilson(k, n)
  sprintf("%-18s %8.4f %8.4f %7.4f %4d/%-3d %5.1f%% [%4.1f, %4.1f]  sd(est) %.4f  mean SE %.4f  SE/sd %.3f",
          lab, truth, mean(est), sd(est) / sqrt(n), k, n, 100 * k / n, 100 * w[1], 100 * w[2],
          sd(est), mean(se), mean(se) / sd(est))
}
cat(sprintf("\n%-18s %8s %8s %7s %9s %6s %s\n", "quantity", "truth", "mean", "mcse", "covered", "rate", "Wilson"))
cat(row("mu intercept", "(Intercept)", 0.4), "\n")
cat(row("mu cond", "cond", 0.9), "\n")
cat(row("log bs", "bs_(Intercept)", log(1.4)), "\n")
cat(row("log sd(mu|s)", "theta_1", log(0.35)), "\n")
cat(row("log sd(lbs|s)", "theta_2", log(0.20)), "\n")

cat("\n-- where the drift intercept's error comes from\n")
est <- get("(Intercept)", "est")
se <- (get("(Intercept)", "upr") - get("(Intercept)", "lwr")) / (2 * qnorm(0.975))
mu_u <- vapply(R, function(r) r$mean_u, 0)
sd_u <- vapply(R, function(r) r$sd_u, 0)
th1 <- get("theta_1", "est")
cat(sprintf("ideal SE of the population mean, 0.35 / sqrt(30): %.4f\n", 0.35 / sqrt(30)))
cat(sprintf("sd of the realized mean of u over replicates: %.4f\n", sd(mu_u)))
cat(sprintf("sd(est) %.4f; sd(est - 0.4 - mean(u)) %.4f (the within-subject part)\n",
            sd(est), sd(est - 0.4 - mu_u)))
cat(sprintf("mean reported SE %.4f; mean exp(theta_1)/sqrt(30) %.4f; mean realized sd(u)/sqrt(30) %.4f\n",
            mean(se), mean(exp(th1)) / sqrt(30), mean(sd_u) / sqrt(30)))
cat(sprintf("cor(est - 0.4, mean(u)) %.3f\n", cor(est - 0.4, mu_u)))
z <- (est - 0.4) / se
cat(sprintf("z = (est - 0.4) / SE: mean %.3f, sd %.3f; |z| > 1.96 on %d of %d\n",
            mean(z), sd(z), sum(abs(z) > 1.96), n))
z2 <- (est - 0.4 - mu_u) / se
cat(sprintf("against the REALIZED mean 0.4 + mean(u): covered on %d of %d\n", sum(abs(z2) < 1.96), n))
k_t <- sum(abs(z) < qt(0.975, 29))
cat(sprintf("with a t(29) critical value %.3f: covered on %d of %d\n", qt(0.975, 29), k_t, n))
cat(sprintf("Shapiro-Wilk on est: p = %.3f\n", shapiro.test(est)$p.value))
saveRDS(data.frame(seed = seeds[!err], est, se, mu_u, sd_u, th1), file.path("dev/phase3b-review2", paste0("cov-summary-", paste(a, collapse = "-"), ".rds")))
