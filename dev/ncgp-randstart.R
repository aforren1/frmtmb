# frm_sample(fit)'s start at the nugget: seed 6601, n 15, true sigma 1e-6 and 1e-3. The
# rule's default against init = "last.par.best". 2 chains in one process, since parallel chains collided on a port at this load.
#   Rscript dev/ncgp-randstart.R <sigma> <mode|default>
a <- commandArgs(TRUE)
s <- as.numeric(a[1])
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(6601)
n <- 15L
x <- round(stats::runif(n, 0, 6), 3)
d <- data.frame(x = x, y = sin(2 * x) + s * stats::rnorm(n))
fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
cat(sprintf("RS fit sigma %.3g warning %s\n", exp(fit$estimates$betad[[1L]]),
            !is.null(gp_nugget_message(fit))))
msg <- character()
args <- list(fit, chains = 2L, cores = 1L, iter = 1000L, warmup = 500L, seed = 11L, refresh = 0,
             control = list(adapt_delta = 0.95, max_treedepth = 10L))
if (a[2] == "mode") args$init <- "last.par.best"
t0 <- proc.time()[[3]]
fs <- withCallingHandlers(suppressWarnings(do.call(frm_sample, args)),
                          message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") })
dr <- posterior::as_draws_df(fs)
cat(sprintf("RS s %.0e %-7s random message %s %.0fs | sigma q50 %.3g rhat %.3f | sdgp q50 %.3g rhat %.3f\n",
            s, a[2], any(grepl("chains start at random", msg)), proc.time()[[3]] - t0,
            stats::median(dr$sigma), posterior::rhat(posterior::extract_variable_matrix(dr, "sigma")),
            stats::median(dr$sdgp_gpx), posterior::rhat(posterior::extract_variable_matrix(dr, "sdgp_gpx"))))
