source("frailty-lib.R")
x <- readLines("frailty-recover.tsv")
rows <- lapply(strsplit(x, "\t"), function(f) {
  kv <- strsplit(f, "=")
  v <- vapply(kv, function(z) trimws(z[2L]), character(1))
  names(v) <- vapply(kv, `[`, character(1), 1L)
  v
})
nm <- names(rows[[which.max(lengths(rows))]])
d <- as.data.frame(do.call(rbind, lapply(rows, function(r) r[nm])),
                   stringsAsFactors = FALSE)
for (j in nm) {
  if (j == "status") next
  d[[j]] <- suppressWarnings(as.numeric(d[[j]]))
}
cat("rows in file", nrow(d), " ok", sum(d$status == "ok"), "\n")
d <- d[d$status == "ok", ]
n <- nrow(d)
f3 <- function(v) formatC(v, digits = 4, format = "g")
truth <- list(beta = 0.6, sd = 0.5, g0 = -1.3 * log(5), g1 = 1.3,
              g2 = 0, g3 = 0)
cat("\n== recovery, ", n, " replicates, mcse at 0.95 = ",
    f3(sqrt(0.95 * 0.05 / n)), " ==\n", sep = "")
cat(sprintf("%-12s %10s %10s %10s %12s\n", "quantity", "truth", "mean",
            "sd", "coverage"))
line <- function(lab, v, tr, cov) {
  cat(sprintf("%-12s %10s %10s %10s %12s\n", lab, f3(tr), f3(mean(v)),
              f3(sd(v)), paste0(sum(cov), "/", n, " ", f3(mean(cov)))))
}
line("beta(trt)", d$beta_frm, truth$beta, d$cov_beta)
line("sd(frailty)", d$sd_frm, truth$sd, d$cov_sd)
line("gamma0", d$g0, truth$g0, d$cov_g0)
line("gamma1", d$g1, truth$g1, d$cov_g1)
line("gamma2", d$g2, truth$g2, d$cov_g2)
line("gamma3", d$g3, truth$g3, d$cov_g3)

ok <- !is.na(d$beta_rst)
cat("\n== agreement with rstpm2 on ", sum(ok), " of ", n,
    " replicates ==\n", sep = "")
bd <- abs(d$beta_frm - d$beta_rst)[ok]
cat("beta:   |diff| mean ", f3(mean(bd)), " max ", f3(max(bd)),
    "; in units of the run's own se: mean ",
    f3(mean(bd / d$se_beta_frm[ok])), " max ",
    f3(max(bd / d$se_beta_frm[ok])), "\n", sep = "")
sr <- abs(d$sd_frm - d$sd_rst)[ok] / d$sd_rst[ok]
cat("sd:     relative |diff| mean ", f3(mean(sr)), " max ", f3(max(sr)),
    "\n", sep = "")
ser <- abs(d$se_beta_frm - d$se_beta_rst)[ok] / d$se_beta_frm[ok]
cat("se(beta): relative |diff| mean ", f3(mean(ser)), " max ",
    f3(max(ser)), "\n", sep = "")
cat("nsx -> RP basis map residual, relative: max ",
    f3(max(d$map_rel[ok])), "\n", sep = "")

cat("\n== the three integration rules ==\n")
lap <- d$ll_frm - d$exact_at_frm
ghq <- (d$ll_rst - d$exact_at_rst)[ok]
gapx <- (d$exact_at_frm - d$exact_at_rst)[ok]
cat("Laplace minus exact, at frmtmb's own optimum: mean ",
    f3(mean(lap)), " range ", f3(min(lap)), " to ", f3(max(lap)),
    "\n", sep = "")
cat("9-node adaptive GH minus exact, at rstpm2's: mean ", f3(mean(ghq)),
    " range ", f3(min(ghq)), " to ", f3(max(ghq)), "\n", sep = "")
cat("exact ll at frmtmb's point minus at rstpm2's: mean ",
    f3(mean(gapx)), " worst ", f3(min(gapx)), "\n", sep = "")
cat("that gap as a fraction of the Laplace offset: mean ",
    f3(mean(abs(gapx) / abs(lap[ok]))), " max ",
    f3(max(abs(gapx) / abs(lap[ok]))), "\n", sep = "")
cat("frmtmb's point is the better one under the EXACT criterion on ",
    sum(gapx > 0), " of ", sum(ok), "\n", sep = "")

cat("\n== convergence and floors ==\n")
cat("conv=0 on ", sum(d$conv == 0), "/", n, ", pdHess on ",
    sum(d$pdHess == 1), "/", n, ", max maxgrad ", f3(max(d$maxgrad)),
    "\n", sep = "")
cat("non-monotone rows: max ", max(d$n_nonmono), " over ", n,
    " fits; deepest censored -log S ", f3(max(d$max_nlogS)),
    " against rp_floored()'s threshold of 19.2\n", sep = "")

cat("\n== the replicates whose optimizer complained ==\n")
bad <- d$conv != 0 | d$maxgrad > 5e-3
cat("flagged ", sum(bad), " of ", n, ": seeds ",
    paste(d$seed[bad], collapse = " "), "\n", sep = "")
if (any(bad)) {
  cat("their agreement with rstpm2, in units of the run's own se: ",
      paste(f3(abs(d$beta_frm - d$beta_rst)[bad] / d$se_beta_frm[bad]),
            collapse = " "), "\n", sep = "")
  cat("their sd relative difference: ",
      paste(f3(abs(d$sd_frm - d$sd_rst)[bad] / d$sd_rst[bad]),
            collapse = " "), "\n", sep = "")
  cat("their maxgrad: ", paste(f3(d$maxgrad[bad]), collapse = " "),
      "\n", sep = "")
}
