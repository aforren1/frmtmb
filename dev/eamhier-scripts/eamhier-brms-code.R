# Lane eamhier: what brms's hierarchical Wiener would be compared
# AGAINST, read off its own generated Stan program rather than argued
# from memory. No compile, no sampling: make_stancode() only.
#
# The check column of item 2.1 asks for a hierarchical Wiener in brms on
# one data set. Before paying for it, this script asks whether the
# comparison is even well posed on the design item 2.1 names, where 20
# of the 30 subjects have a true non-decision time above the GLOBAL
# fastest response.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-brms-code.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("brms", format(packageVersion("brms")), "| rstan",
    format(packageVersion("rstan")), "| StanHeaders",
    format(packageVersion("StanHeaders")), "\n")

d <- eamhier_data(20260908L, "A")
t0 <- attr(d, "ndt_subject")
floors <- as.numeric(tapply(d$rt, d$s, min))
cat("global min(rt)          ", min(d$rt), "\n")
cat("true ndt range          ", min(t0), max(t0), "\n")
cat("subjects with true ndt above the GLOBAL bound:",
    sum(t0 > min(d$rt)), "of", length(t0), "\n")
cat("subjects with true ndt above their OWN bound: ",
    sum(t0 > floors), "of", length(t0), "\n")

code <- brms::make_stancode(
  brms::bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
           ndt ~ 1 + (1 | s), bias = 0.5),
  family = brms::wiener(), data = d)
cat("\n== the ndt lines of brms's generated program ==\n")
ln <- strsplit(as.character(code), "\n", fixed = TRUE)[[1L]]
keep <- grep("ndt|min_Y|Y\\[", ln)
cat(paste(ln[keep], collapse = "\n"), "\n")

sd <- brms::make_standata(
  brms::bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
           ndt ~ 1 + (1 | s), bias = 0.5),
  family = brms::wiener(), data = d)
cat("\nmin_Y that brms hands Stan:", sd$min_Y, "\n")
