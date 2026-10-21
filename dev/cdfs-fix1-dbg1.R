# Fix round 1 debugging: the review's grad5 interval rows (skew_normal,
# alpha = -10) still show AD != FD. Per row: the interval's log mass by
# the new log_window(), by sn::psn() at high precision (Rmpfr-free: the
# integral of the density), and AD vs FD of each row in each parameter.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
src <- readLines("dev/cdfs-rev-grad.R")
cut <- grep("^pick <- ", src)
Sys.setenv(GROUP = 2)
eval(parse(text = src[1:(cut - 1)]))
fm <- fams[["skew_normal"]]
dpn <- setdiff(frmtmb:::as_frmtmb_family(fm$fam)$dpars, "mu")
fo <- do.call(bf, c(list(y | cens(cc, y2) ~ x),
                    lapply(dpn, function(dp) as.formula(paste(dp, "~ 1")))))
d <- fm$d
d$cc[d$cc != 2] <- 0
ob <- frm(fo, family = fm$fam, data = d, dry_run = "objective")$obj
p <- ob$par
jd <- which(names(p) == "betad")
p[jd[match("alpha", dpn)]] <- -10
g <- as.numeric(ob$gr(p))
fd <- numDeriv::grad(ob$fn, p, method = "Richardson")
print(rbind(names(p), signif(g, 8), signif(fd, 8)))
# FD at fixed steps
for (h in c(1e-3, 1e-4, 1e-5, 1e-6)) {
  fdh <- vapply(seq_along(p), function(j) {
    e <- replace(numeric(length(p)), j, h)
    (ob$fn(p + e) - ob$fn(p - e)) / (2 * h)
  }, 0)
  cat("h", h, ":", signif(fdh, 8), "\n")
}
cat("END\n")
