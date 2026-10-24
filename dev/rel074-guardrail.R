# Guardrail 1 on the merged build (dev/core-boundary.md): every brms
# 2.23.0 name, probed with the probes and the classifier of
# tests/testthat/test-guardrail1.R, first in a child with core alone,
# then, for each extension a refusal names, in a child with that
# extension loaded. Prints every name that is not "core" with its class
# in both arms, so the record shows cox() resolving in frmtmb.survival
# and car() and spde() in frmtmb.spatial. Adapted from
# dev/spatial-guard-report.R.
#   Rscript dev/rel074-guardrail.R <lib;lib>
libs <- strsplit(commandArgs(trailingOnly = TRUE)[1], ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
src <- readLines("tests/testthat/test-guardrail1.R")
first <- grep("^test_that", src)[1]
eval(parse(text = paste(src[seq_len(first - 1)], collapse = "\n")))
pr <- guard_probes()
res <- guard_run(pr)
stopifnot(isTRUE(attr(res, "ok")))
cls <- mapply(guard_classify, res, pr$name)
cat("probes:", nrow(pr), "\n")
print(table(pr$kind, cls))
homes <- c(frmtmb:::extension_family_homes, frmtmb:::extension_special_homes,
           frmtmb:::extension_aterm_homes)
cat("home tables:", paste(names(homes), homes, sep = " = ", collapse = "; "),
    "\n")
cat("\nnames that are not core, core alone:\n")
nc <- which(cls != "core")
for (i in nc) {
  cat(sprintf("  %-8s %-34s %-22s %s\n", pr$kind[i], pr$name[i], cls[i],
              substr(res[i], 1, 100)))
}
cat("\nthe extension arm:\n")
for (pkg in unique(sub("^ext:", "", cls[startsWith(cls, "ext:")]))) {
  at <- which(cls == paste0("ext:", pkg))
  if (!nzchar(system.file(package = pkg))) {
    cat("  ", pkg, ": NOT INSTALLED\n")
    next
  }
  r2 <- guard_run(pr[at, , drop = FALSE], pkg)
  c2 <- mapply(guard_classify, r2, pr$name[at])
  for (j in seq_along(at)) {
    cat(sprintf("  with %-16s %-8s %-34s %-8s %s\n", pkg, pr$kind[at[j]],
                pr$name[at[j]], c2[j], substr(r2[j], 1, 60)))
  }
}
# spde() is no brms special, so the vocabulary has no probe for it; its
# refusal by name is probed from car()'s plain spelling, both arms
sp <- pr[pr$name == "car" & pr$ext == "ok", , drop = FALSE][1, ]
sp$name <- "spde"
sp$code <- gsub("car(", "spde(", sp$code, fixed = TRUE)
cat("\nspde probe:", sp$code, "\n")
for (pk in list(character(), "frmtmb.spatial")) {
  r3 <- guard_run(sp, pk)
  cat(sprintf("  with %-16s spde: %s\n", if (length(pk)) pk else "core alone",
              substr(r3[1], 1, 110)))
}
# the bhaz() refusal on the merged build
cat("\nbhaz():", tryCatch(frmtmb:::check_unsupported_aterm("bhaz"),
                         error = conditionMessage), "\n")
