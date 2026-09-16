## The S3 generic CONTRACT scan, and the table it produces.
##
## The lane's rule was "refuse what the METHOD does not have". R's S3
## contract is "a method tolerates the arguments its GENERIC documents",
## which is why every stats default method carries `use.fallback`,
## `env`, `data`, `na.rm` or `trace` and ignores them. Three rounds of
## `.allow` patches each uncovered the next hop; this finds them all at
## once, and `tests/testthat/test-arg-refusal.R` runs the same scan so
## the NEXT one is found by the suite.
##
## For every S3 generic either package registers a dots-refusing method
## on, walk the ASTs of the base packages and the interop packages, find
## real CALLS to that generic, and report the named arguments our method
## has no formal for. Those are exactly the calls the refusal breaks.
##
## Run: Rscript dev/argspell-contract.R
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

source(file.path("dev", "argspell-contract-scan.R"))

pkgs <- c("frmtmb", "frmtmb.sample")
# R's OWN packages only. The S3 contract is R's; a third-party
# package's calls to `plot()` on its own objects are not a claim about
# what our method must tolerate, and scanning them buried the signal.
scan_pkgs <- c("stats", "base", "utils", "graphics", "grDevices",
               "methods")
scan_pkgs <- Filter(function(p) requireNamespace(p, quietly = TRUE),
                    scan_pkgs)

# With the shipped table, and then without it, so the run shows both
# what remains and what the table is carrying.
acc <- frm_contract_scan(pkgs, scan_pkgs,
                         contract = frmtmb:::s3_contract_args)
bare <- frm_contract_scan(pkgs, scan_pkgs, contract = NULL)
cat("---- GENERATED: S3 generic contract scan ----\n")
cat("packages scanned:", paste(scan_pkgs, collapse = ", "), "\n")
cat(sprintf("generics with a dots-refusing method: %d\n", acc$n_generics))
cat(sprintf("call sites, table EMPTY             : %d\n",
            nrow(bare$hits)))
cat(sprintf("call sites, table APPLIED           : %d\n",
            nrow(acc$hits)))
h0 <- bare$hits[order(bare$hits$generic, bare$hits$arg), ]
cat("\nwhat the table carries, generic to the names R itself passes:\n")
ag0 <- tapply(h0$arg, h0$generic,
              function(z) paste(sort(unique(z)), collapse = ", "))
for (g in sort(names(ag0))) cat(sprintf("  %-14s %s\n", g, ag0[[g]]))
if (nrow(acc$hits)) {
  cat("\nSTILL UNCOVERED:\n")
  h <- acc$hits[order(acc$hits$generic, acc$hits$arg), ]
  for (i in seq_len(nrow(h))) {
    cat(sprintf("  %-14s %-22s %s\n", h$generic[i], h$arg[i],
                h$where[i]))
  }
}
cat("---- END GENERATED: S3 generic contract scan ----\n")
