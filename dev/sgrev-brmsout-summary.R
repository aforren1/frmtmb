root <- "dev/sgrev-out/brmsout"
ctl1 <- readRDS(file.path(root, "none1.rds"))
ctl2 <- readRDS(file.path(root, "none2.rds"))
nms <- names(ctl1)

cat("== reviewer output harness, dev/sgrev-brmsout.R ==\n")
cat("cached 400-row lognormal brmsfit, seed 20260915 before every",
    "call\n")
cat("60 calls. Each compared on VALUE, WARNINGS, MESSAGES and the",
    "PRINTED form.\n\n")

# determinism control: two brms-only processes, different library
# paths but the same brms
dif <- function(a, b, fields = c("val", "warn", "msg", "printed")) {
  out <- character()
  for (nm in names(a)) {
    x <- a[[nm]]; y <- b[[nm]]
    bad <- fields[vapply(fields, function(f)
      !identical(x[[f]], y[[f]]), NA)]
    if (length(bad)) out <- c(out, paste0(nm, "(", paste(bad,
      collapse = ","), ")"))
  }
  out
}
d0 <- dif(ctl1, ctl2)
cat("CONTROL brms-only vs brms-only, two processes:",
    length(d0), "of", length(nms), "differ\n")
if (length(d0)) cat("   ", paste(d0, collapse = " "), "\n")
cat("   (these calls are not deterministic across processes and are",
    "excluded below)\n\n")
stable <- setdiff(nms, sub("[(].*$", "", d0))
cat("deterministic calls:", length(stable), "of", length(nms), "\n\n")

arms <- c("BASE-S", "BASE-T", "FIX-S", "FIX-T", "FIX-U", "FIX-SC")
cat(sprintf("%-8s %-28s %-9s %-9s %-9s %-9s\n", "arm", "",
            "value", "warnings", "messages", "printed"))
for (a in arms) {
  f <- file.path(root, paste0(a, ".rds"))
  if (!file.exists(f)) next
  x <- readRDS(f)
  cnt <- function(fld) sum(vapply(stable, function(nm)
    !identical(x[[nm]][[fld]], ctl1[[nm]][[fld]]), NA))
  cat(sprintf("%-8s %-28s %-9s %-9s %-9s %-9s\n", a, "",
              paste0(cnt("val"), "/", length(stable)),
              paste0(cnt("warn"), "/", length(stable)),
              paste0(cnt("msg"), "/", length(stable)),
              paste0(cnt("printed"), "/", length(stable))))
  d <- dif(x[stable], ctl1[stable])
  if (length(d)) {
    cat("   differing: ", paste(d, collapse = " "), "\n")
    for (nm in sub("[(].*$", "", d)) {
      if (!identical(x[[nm]]$warn_txt, ctl1[[nm]]$warn_txt))
        cat("     ", nm, " warnings control=[",
            paste(ctl1[[nm]]$warn_txt, collapse = " | "), "] arm=[",
            paste(x[[nm]]$warn_txt, collapse = " | "), "]\n", sep = "")
      if (!identical(x[[nm]]$msg_txt, ctl1[[nm]]$msg_txt))
        cat("     ", nm, " messages control=[",
            paste(ctl1[[nm]]$msg_txt, collapse = " | "), "] arm=[",
            paste(x[[nm]]$msg_txt, collapse = " | "), "]\n", sep = "")
      if (!identical(x[[nm]]$err_txt, ctl1[[nm]]$err_txt))
        cat("     ", nm, " error control=[", ctl1[[nm]]$err_txt,
            "] arm=[", x[[nm]]$err_txt, "]\n", sep = "")
      if (!identical(x[[nm]]$cls, ctl1[[nm]]$cls))
        cat("     ", nm, " class control=", ctl1[[nm]]$cls,
            " arm=", x[[nm]]$cls, "\n", sep = "")
    }
  }
}
cat("\n-- every call's ERROR text in the brms-only control --\n")
for (nm in nms) if (identical(ctl1[[nm]]$kind, "ERROR"))
  cat(sprintf("  %-22s %s\n", nm, substr(ctl1[[nm]]$err_txt, 1, 90)))
cat("\n-- calls that WARN in the brms-only control --\n")
for (nm in nms) if (length(ctl1[[nm]]$warn_txt))
  cat(sprintf("  %-22s %s\n", nm,
              substr(paste(ctl1[[nm]]$warn_txt, collapse = " | "), 1, 90)))
cat("\n-- calls that MESSAGE in the brms-only control --\n")
for (nm in nms) if (length(ctl1[[nm]]$msg_txt))
  cat(sprintf("  %-22s %s\n", nm,
              substr(paste(ctl1[[nm]]$msg_txt, collapse = " | "), 1, 90)))
cat("DONE\n")
