## frm_sample() of y ~ mo(k), k in 1..4, against brms 2.23.0 on the same
## data (seed 5, n = 300); both 2 chains of 2000 on seed 3. The numbers
## the gated test in frmtmb.sample's test-mo-coding.R asserts.
##   Rscript dev/mocode-agree.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-mocode/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("== frmtmb", format(packageVersion("frmtmb")), arm,
    dirname(system.file(package = "frmtmb")), "\n")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
code <- brms::stancode(y ~ mo(k), data = d)
sdat <- brms::standata(y ~ mo(k), data = d)
f <- tempfile(fileext = ".stan")
writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
path <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"),
                  paste0(unname(tools::md5sum(f)), ".rds"))
mod <- if (file.exists(path)) readRDS(path) else {
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
sf <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 2,
                                       iter = 2000, seed = 3, refresh = 0))
mb <- as.matrix(sf, pars = c("b_Intercept", "bsp", "sigma", "simo_1"))
colnames(mb) <- c("b_Intercept", "bsp_mok", "sigma",
                  paste0("simo_mok1[", 1:3, "]"))
w <- character(0)
ds <- withCallingHandlers(suppressMessages(frm_sample(
  y ~ mo(k), data = d, chains = 2, iter = 2000, seed = 3, refresh = 0)),
  warning = function(w0) {
    w <<- c(w, conditionMessage(w0))
    invokeRestart("muffleWarning")
  })
cat("frm_sample warnings:", length(w), "\n")
if (length(w)) cat(paste0("  ", substr(w, 1, 100)), sep = "\n")
mf <- as.matrix(ds)
cat("frmtmb simo variables:", grep("^simo", colnames(mf), value = TRUE),
    "\n")
mc <- function(v) posterior::mcse_mean(matrix(v, ncol = 2))
cat(sprintf("%-14s %9s %9s %9s %9s %8s %8s\n", "par", "mean_frm",
            "mean_brms", "sd_frm", "sd_brms", "z", "sd_ratio"))
for (p in colnames(mb)) {
  z <- (mean(mf[, p]) - mean(mb[, p])) / sqrt(mc(mf[, p])^2 + mc(mb[, p])^2)
  cat(sprintf("%-14s %9.4f %9.4f %9.4f %9.4f %8.2f %8.3f\n", p,
              mean(mf[, p]), mean(mb[, p]), sd(mf[, p]), sd(mb[, p]), z,
              sd(mf[, p]) / sd(mb[, p])))
}
cat("brms rhat max:", max(apply(mb, 2, function(v)
  posterior::rhat(matrix(v, ncol = 2)))), "\n")
