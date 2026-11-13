# Lane nlstop: one row of the review's fresh sets in detail. Fits the
# default, prints its estimates, warnings and messages, the SE check's
# verdict, the best of the 16-17 reference starts and where it ends.
# Usage: NLSTOP_ARM=before|after Rscript dev/nlstop-row.R <set> <fam> <body> <seed>
#   set: fresh (dev/nlstop-fresh.R) or fresh2 (dev/nlstop-fresh2.R)
a <- commandArgs(TRUE)
set <- a[1]; fam <- a[2]; body <- a[3]; seed <- as.numeric(a[4])
source("dev/nlstop-common.R")
src <- parse(sprintf("dev/nlstop-%s.R", set))
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
d <- gen(fam, body, seed * 1000 + match(body, names(B)))
fm <- form(fam, body)
show <- function(lab, r) {
  if (iserr(r)) { cat(lab, "ERROR", r$msg, "\n"); return(invisible()) }
  cat(sprintf("%s: ll %.6f llp %.6f code %s (%s)\n", lab, llv(r), llp(r),
              r$opt$convergence, r$opt$message))
  b <- r$estimates$beta
  cat("  beta:", paste(sprintf("%s=%.6g", names(b), b), collapse = "; "), "\n")
  for (i in seq_along(attr(r, "w"))) {
    cat("  W[", attr(r, "wc")[i], "]:",
        substr(gsub("\n", " ", attr(r, "w")[i]), 1, 300), "\n")
  }
  for (i in seq_along(attr(r, "m"))) {
    cat("  M[", attr(r, "mc")[i], "]:",
        substr(gsub("\n", " ", attr(r, "m")[i]), 1, 300), "\n")
  }
}
r <- tt(frm(fm, family = F[[fam]], data = d))
show("default", r)
ro <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                          tt(frm(fm, family = F[[fam]], data = d)))
show("probe off", ro)
nm <- names(r$estimates$beta)
np <- sub("_.*$", "", nm)
sc <- B[[body]]$sc[np]
sc[is.na(sc)] <- 0.5
set.seed(seed * 13 + 7)
starts <- lapply(1:16, function(i) {
  setNames(rnorm(length(nm)) * sc * if (i > 8) 3 else 1, nm)
})
if (!is.null(B[[body]]$tr)) {
  tr <- setNames(rep(0, length(nm)), nm)
  ic <- grepl("_\\(Intercept\\)$", nm)
  tr[ic] <- B[[body]]$tr[np[ic]]
  starts <- c(starts, list(tr))
}
vals <- vapply(starts, function(st) {
  llv(tt(frm(fm, family = F[[fam]], data = d, se = FALSE,
             start = list(beta = st))))
}, 0)
cat("reference starts:", format(sort(vals, decreasing = TRUE), digits = 9),
    "\n")
ib <- which.max(vals)
rb <- tt(frm(fm, family = F[[fam]], data = d, start = list(beta = starts[[ib]])))
show(sprintf("best start %d", ib), rb)
cat("  start:", format(starts[[ib]], digits = 4), "\n")
