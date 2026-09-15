# Lane eamhier: what the post-fit surface calls things, at a size that
# costs nothing. It exists so that the replicate script reads
# coefficients by name rather than by position.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-smoke.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

d <- eamhier_data(20260910L, "C", ns = 6L, nt = 80L)
cat("rows", nrow(d), "min rt", min(d$rt), "\n")

t0 <- Sys.time()
fit <- frm(eamhier_form("C", "pg"), family = eamhier_family("C"),
           data = d, se = TRUE)
cat("fit_s", as.numeric(difftime(Sys.time(), t0, units = "secs")), "\n")

cat("== fixef ==\n"); print(unlist(fixef(fit)))
cat("== confint rownames ==\n")
ci <- suppressWarnings(stats::confint(fit))
print(ci)
cat("== VarCorr ==\n"); print(VarCorr(fit))
cat("== opt$par names ==\n"); print(names(fit$opt$par))
cat("== sdr ==\n"); print(summary(fit$sdr, "fixed"))
cat("== diagnose ==\n")
dg <- frmtmb::diagnose(fit, quiet = TRUE)
print(names(dg))
str(dg[c("convergence", "max_grad", "pdHess")])
cat("== ndt_time ==\n")
key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
print(as.numeric(ndt_time(fit, newdata = key)))
print(attr(d, "ndt_subject"))
