# Part B: fixture 1 refit with count * (1 + 2e-15 * j), j = -10..10, so
# that the stop moves by rounding; the verdicts the new test in
# test-se-conv.R asserts, per draw.
#   CODE8_ARM=before|after Rscript dev/code8-fix1-draws.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/code8-common.R")
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
c0 <- d$count
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
rows <- list()
for (j in -10:10) {
  d$count <- c0 * (1 + 2e-15 * j)
  r <- rfit(form, family = student(), data = d)
  f <- r$f
  nm <- fns$outer_par_names(f)
  sdr <- quiet(fns$sdr_of(f))
  se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  fi <- tryCatch(quiet(fitted(f, dpar = "sigma")), error = function(e) NULL)
  rows[[length(rows) + 1L]] <- data.frame(
    j = j, obj = f$opt$objective, code = f$opt$convergence,
    theta3 = f$opt$par[match("theta_3", nm)],
    se_sig = paste(format(se[match(c("sigma_(Intercept)", "sigma_Trt1"), nm)],
                          digits = 4), collapse = ","),
    h326 = !is.null(fi) && isTRUE(all(fi > 0)),
    grad_w = any(grepl("Large maximum absolute gradient", r$w)),
    short_w = any(grepl("short of the maximum|curvature at the stop", r$w)),
    bnd_m = sum(grepl("frmtmb_boundary_fit", r$mc)),
    lost = paste(sprintf("%s:%s", names(sdr$se_lost), sdr$se_lost),
                 collapse = " "))
}
X <- do.call(rbind, rows)
write.table(X, sprintf("dev/code8-out/fix1-draws-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
X$obj <- sprintf("%.7f", X$obj)
print(X, right = FALSE)
cat("h326 TRUE", sum(X$h326), "of", nrow(X), "; grad warning", sum(X$grad_w),
    "; short warning", sum(X$short_w), "; boundary message",
    sum(X$bnd_m > 0), "\n")
cat("DONE\n")
