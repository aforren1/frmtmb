# Lane eamhier: check_laplace() on one data set, which is the half of
# item 2.1's check column that needs no Stan program of its own.
#
# It samples the density the fit maximized with tmbstan, so it compares
# NUTS against the Laplace mode and the Wald standard errors on the SAME
# likelihood. That is the question a hierarchical Wald interval raises:
# the coverage this lane measures is a Wald coverage, and if it misses,
# the Laplace approximation is the first suspect.
#
# It runs the joint density with the random effects sampled, so its cost
# is set by the number of leapfrog steps and not by the fit's, which is
# why the size is an argument and why the pilot is small.
#
# Run:
#   Rscript --vanilla dev/eamhier-scripts/eamhier-laplace.R \
#     <arm A|C> <seed> <ns> <nt> <chains> <iter> <outfile>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L)
arm <- args[[1L]]
seed <- as.integer(args[[2L]])
ns <- as.integer(args[[3L]])
nt <- as.integer(args[[4L]])
chains <- as.integer(args[[5L]])
iter <- as.integer(args[[6L]])
out <- args[[7L]]

source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.sample)
})
cat("StanHeaders", format(packageVersion("StanHeaders")),
    "| rstan", format(packageVersion("rstan")),
    "| tmbstan", format(packageVersion("tmbstan")), "\n")
# The silent failure the lane rules name: a tmbstan built against
# StanHeaders 2.39 samples a standard normal and says nothing.
cat("tmbstan build broken:",
    frmtmb.sample:::tmbstan_build_broken(), "\n")

d <- eamhier_data(seed, arm, ns = ns, nt = nt)
t0 <- Sys.time()
fit <- frm(eamhier_form(arm, "pg"), family = eamhier_family(arm),
           data = d, se = TRUE)
cat("fit_s", as.numeric(difftime(Sys.time(), t0, units = "secs")),
    "logLik", as.numeric(logLik(fit)), "\n")

t1 <- Sys.time()
cl <- check_laplace(fit, chains = chains, iter = iter, refresh = 0)
cl_s <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
cat("check_laplace_s", cl_s, "\n")
cl$arm <- arm
cl$seed <- seed
cl$ns <- ns
cl$nt <- nt
cl$chains <- chains
cl$iter <- iter
cl$check_laplace_s <- cl_s
print(cl, digits = 4)
utils::write.table(cl, out, sep = "\t", row.names = FALSE,
                   quote = FALSE)
