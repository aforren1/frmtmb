root <- "dev/sgrev-out/drawsout"
ctl <- readRDS(file.path(root, "BASE-N.rds"))
ctl2 <- readRDS(file.path(root, "BASE-N2.rds"))
nms <- names(ctl)
flds <- c("val", "warn", "msg", "printed")
dif <- function(a, b) {
  out <- character()
  for (nm in names(a)) {
    bad <- flds[vapply(flds, function(f)
      !identical(a[[nm]][[f]], b[[nm]][[f]]), NA)]
    if (length(bad)) out <- c(out, paste0(nm, "(",
      paste(bad, collapse = ","), ")"))
  }
  out
}
cat("== draws-side output, BASE against FIX, dev/sgrev-drawsout.R ==\n")
cat("46 calls on one cached frmtmb_draws (4 chains x 500), seed",
    "20260915 before each\n")
cat("reference arm: BASE, frmtmb.sample alone\n\n")
d0 <- dif(ctl2, ctl)
cat("CONTROL, the reference arm run twice: ", length(d0), " of ",
    length(nms), " differ\n", sep = "")
if (length(d0)) cat("   ", paste(d0, collapse = " "), "\n")
stable <- setdiff(nms, sub("[(].*$", "", d0))
cat("deterministic calls: ", length(stable), " of ", length(nms),
    "\n\n", sep = "")

for (a in c("BASE-S", "BASE-T", "BASE-P", "BASE-Q",
            "FIX-N", "FIX-S", "FIX-T", "FIX-P", "FIX-Q")) {
  f <- file.path(root, paste0(a, ".rds"))
  if (!file.exists(f)) next
  x <- readRDS(f)
  d <- dif(x[stable], ctl[stable])
  cat(sprintf("%-8s differs on %2d of %d calls\n", a, length(d),
              length(stable)))
  for (e in d) {
    nm <- sub("[(].*$", "", e)
    cat("   ", e, "\n")
    if (!identical(x[[nm]]$warn_txt, ctl[[nm]]$warn_txt))
      cat("      warn base=[", paste(ctl[[nm]]$warn_txt, collapse = " | "),
          "] arm=[", paste(x[[nm]]$warn_txt, collapse = " | "), "]\n",
          sep = "")
    if (!identical(x[[nm]]$msg_txt, ctl[[nm]]$msg_txt))
      cat("      msg  base=[", paste(ctl[[nm]]$msg_txt, collapse = " | "),
          "] arm=[", paste(x[[nm]]$msg_txt, collapse = " | "), "]\n",
          sep = "")
    if (!identical(x[[nm]]$err_txt, ctl[[nm]]$err_txt))
      cat("      err  base=[", substr(ctl[[nm]]$err_txt, 1, 100),
          "] arm=[", substr(x[[nm]]$err_txt, 1, 100), "]\n", sep = "")
  }
}
cat("\n-- the reference arm's own errors, for context --\n")
for (nm in nms) if (identical(ctl[[nm]]$kind, "ERROR"))
  cat(sprintf("  %-22s %s\n", nm, substr(ctl[[nm]]$err_txt, 1, 95)))
cat("DONE\n")
