# What keeping the refit data on the fit costs, and what recovering it
# from the workspace (the call's `data` symbol) gives when that object
# changes. Arg: base | lane. Log: dev/drop1raw-out/memory-<arm>.txt
# Data: seed 3, n = 1e5, 50 unused numeric columns beside the 5 the
# model reads.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(3)
n <- 1e5
d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1, x3 = runif(n),
                g = factor(sample(sprintf("g%03d", 1:200), n, TRUE)))
d$z <- 1 + d$x + log(d$x2) + sin(3 * d$x3) + rnorm(200)[d$g] + rnorm(n)
for (j in 1:50) d[[sprintf("u%02d", j)]] <- rnorm(n)
mb <- function(x) as.numeric(object.size(x)) / 2^20
ser <- function(x, gz) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf))
  saveRDS(x, tf, compress = gz)
  file.size(tf) / 2^20
}
fit <- q(frm(bf(z ~ x + log(x2) + splines::ns(x3, 3) + (1 | g)), data = d))
fr <- fit$frame
cat(sprintf("rows %d, data columns %d\n", nrow(d), ncol(d)))
cat(sprintf("%-34s %10s %12s %12s\n", "object", "in memory", "rds raw",
            "rds gzip"))
row <- function(lab, x) {
  cat(sprintf("%-34s %8.2f MB %9.2f MB %9.2f MB\n", lab, mb(x),
              ser(x, FALSE), ser(x, TRUE)))
}
row("the data, all 55 columns", d)
row("the data, the 5 model columns", d[c("z", "x", "x2", "x3", "g")])
row("fit$frame$data_frame", fr$data_frame)
row("fit$frame$raw_vars", fr$raw_vars)
cat("raw_vars columns:", paste(names(fr$raw_vars), collapse = ","), "\n")
fit_light <- fit
fit_light$frame$raw_vars <- NULL
row("whole fit", fit)
row("whole fit without raw_vars", fit_light)
if (exists("fit_refit_data", asNamespace("frmtmb"))) {
  rd <- frmtmb:::fit_refit_data(fit)
  row("fit_refit_data(fit), built on call", rd)
  cat("refit data columns:", paste(names(rd), collapse = ","), "\n")
  # a vector the model frame holds is not copied into the refit data
  cat("refit data shares the model frame's z:",
      identical(.Internal(address(rd$z)),
                .Internal(address(fr$data_frame$z))), "\n")
}

cat("\n== the data object changes after the fit (n = 2000 subset) ==\n")
d <- d[1:2000, c("z", "x", "x2", "x3", "g")]
fit <- q(frm(bf(z ~ x + log(x2)), data = d))
ll0 <- as.numeric(logLik(fit))
d$x2 <- d$x2^2 + 1
u <- q(update(fit))
cat(sprintf("x2 replaced by x2^2 + 1 in d: update(fit) logLik %.10f, fit %.10f, same %s\n",
            as.numeric(logLik(u)), ll0,
            identical(as.numeric(logLik(u)), ll0)))
d <- d[-(1:10), ]
u <- q(update(fit))
cat(sprintf("10 rows deleted from d: update(fit) nobs %d, fit nobs %d\n",
            nobs(u), nobs(fit)))
