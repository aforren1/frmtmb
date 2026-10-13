# Punch round 1, m3: test-reparam.R's correlated-model agreement gate,
# measured in Monte Carlo standard errors over sampler seeds. Per seed
# and iteration count, the z of every compared mean and sd:
#   z = (est_centered - est_noncentered) / sqrt(mcse_c^2 + mcse_n^2)
# with each chain's own MCSE (posterior::mcse_mean, posterior::mcse_sd).
#   Rscript dev/moparity-rp-gate.R base|lane first last iter[,iter]
# Output: dev/moparity-log/rp-gate-<arm>-<first>-<iters>.rds
args <- commandArgs(TRUE)
arm <- args[1]
seeds <- as.integer(args[2]):as.integer(args[3])
iters <- as.integer(strsplit(args[4], ",")[[1]])
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"),
    "frmtmb.sample from", find.package("frmtmb.sample"), "\n")
# the data design: RP_DESIGN = "ng,n,sd1,sd2"; the default is the
# test's original data (12 groups, 120 rows, sds 0.7 and 0.4)
des <- as.numeric(strsplit(Sys.getenv("RP_DESIGN", "12,120,0.7,0.4"),
                           ",")[[1]])
tag <- gsub(",", "_", Sys.getenv("RP_DESIGN", "orig"))
cat("design", des, "\n")
set.seed(31)
ng <- as.integer(des[1])
nn <- as.integer(des[2])
dd <- data.frame(x = stats::rnorm(nn),
                 g = factor(rep(seq_len(ng), length.out = nn)))
u <- matrix(stats::rnorm(2 * ng), 2)
# a fifth value is the correlation of the two random effects
if (length(des) >= 5L) u[2, ] <- des[5] * u[1, ] + sqrt(1 - des[5]^2) * u[2, ]
dd$y <- stats::rnorm(nn, 1 + 0.5 * dd$x + des[3] * u[1, dd$g] +
                       des[4] * u[2, dd$g] * dd$x, 1)
zs <- function(cen, ncp) {
  nms <- setdiff(colnames(cen), "lp__")
  out <- list()
  for (nm in nms) {
    a <- cen[, nm]
    b <- ncp[, nm]
    zm <- (mean(a) - mean(b)) /
      sqrt(posterior::mcse_mean(a)^2 + posterior::mcse_mean(b)^2)
    zsd <- (stats::sd(a) - stats::sd(b)) /
      sqrt(posterior::mcse_sd(a)^2 + posterior::mcse_sd(b)^2)
    out[[nm]] <- data.frame(var = nm, z_mean = zm, z_sd = zsd,
                            ess_c = posterior::ess_mean(a),
                            ess_n = posterior::ess_mean(b),
                            sd_c = stats::sd(a), sd_n = stats::sd(b))
  }
  do.call(rbind, out)
}
res <- list()
for (s in seeds) {
  for (it in iters) {
    mk <- function(ncp) {
      suppressWarnings(suppressMessages(
        frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
                   chains = 1, iter = it, refresh = 0, seed = s,
                   reparameterize = ncp)))
    }
    cen <- mk(FALSE)
    ncp <- mk(TRUE)
    z <- zs(cen$draws, ncp$draws)
    z$seed <- s
    z$iter <- it
    res[[length(res) + 1L]] <- z
    hyp <- z[!grepl("^r_", z$var), ]
    cat(sprintf(paste("seed %2d iter %4d: max|z| mean %.2f sd %.2f",
                      "(hyper %.2f %.2f) min ess centered %.0f\n"),
                s, it, max(abs(z$z_mean)), max(abs(z$z_sd)),
                max(abs(hyp$z_mean)), max(abs(hyp$z_sd)), min(z$ess_c)))
  }
}
out <- do.call(rbind, res)
saveRDS(out, sprintf("dev/moparity-log/rp-gate-%s-%s-%d-%s.rds", tag, arm,
                     seeds[1], paste(iters, collapse = "_")))
