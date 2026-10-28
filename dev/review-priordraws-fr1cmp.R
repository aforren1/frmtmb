# Reviewer: mi(sdy) "only" brms vs frmtmb (fr1-mi-*.rds), z-scores of
# mean and quantiles with MCSE; the latent Yl[1] against normal(y[1], 0.4).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
L <- "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-log"
b <- readRDS(file.path(L, "fr1-mi-brms.rds")); f <- readRDS(file.path(L, "fr1-mi-frm.rds"))
st <- function(x) { m <- matrix(x, ncol = 4)
  c(mean(x), posterior::mcse_mean(m), sapply(c(.05, .5, .95), function(p)
    c(unname(quantile(x, p)), unname(posterior::mcse_quantile(m, p))))) }
zs <- function(lab, xb, xf) { a <- st(xb); c2 <- st(xf)
  z <- (a[c(1, 3, 5, 7)] - c2[c(1, 3, 5, 7)]) / sqrt(a[c(2, 4, 6, 8)]^2 + c2[c(2, 4, 6, 8)]^2)
  cat(sprintf("%-12s brms mean %7.3f | frm mean %7.3f | z %s\n", lab, a[1], c2[1],
              paste(sprintf("%5.2f", z), collapse = " "))); z }
z <- c()
for (v in c("b_Intercept", "b_x", "sigma")) z <- c(z, zs(v, b$draws[[v]], f$draws[[v]]))
if ("Yl[1]" %in% names(b$draws)) z <- c(z, zs("Yl[1]/miss_1", b$draws[["Yl[1]"]], f$yl1))
set.seed(1); x <- rnorm(40); y1 <- 1 + x[1] + rnorm(40)[1]
for (j in 1:3) z <- c(z, zs(paste0("yrep[,", j, "]"), b$yrep[, j], f$yrep[, j]))
cat("frm latent 1: mean", mean(f$yl1), "sd", sd(f$yl1), "(exact sd 0.4)\n")
cat(sprintf("SUMMARY %d z, max |z| %.2f\n", length(z), max(abs(z))))
