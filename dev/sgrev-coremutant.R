# Item 7, first half. The lane's new core block fails on BASE with
# "unused argument", which is the WEAK form: the symbol is absent, not
# the behavior wrong. Construct the BEHAVIOURAL failure: a core that
# TAKES the argument and does not validate it. If the block still
# fails there, it is testing the refusal and not the signature.
MUT <- "C:/Users/adf44/source/r/sgrev-mut-novalidate"
MUTLIB <- "C:/Users/adf44/source/r/sgrev-mut-novalidate-lib"
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen"
unlink(MUT, recursive = TRUE)
dir.create(MUT, recursive = TRUE)
for (p in c("DESCRIPTION", "NAMESPACE", "R", "src", "man", "inst",
            "tests", "data")) {
  from <- file.path(src, p)
  if (file.exists(from)) file.copy(from, MUT, recursive = TRUE)
}
f <- file.path(MUT, "R", "generic-owners.R")
s <- readChar(f, file.size(f), useBytes = TRUE)
old <- paste0('  if (!ok) {\n',
              '    stop("frm_install_generics(owners =) must be a named list of ",\n',
              '         "non-empty character vectors: one entry per generic name, ",\n',
              '         "naming the packages that own it in the order to prefer",\n',
              '         call. = FALSE)\n  }')
n <- length(gregexpr(old, s, fixed = TRUE)[[1]])
cat("the validation block was found exactly once: ",
    n == 1L && gregexpr(old, s, fixed = TRUE)[[1]][1] > 0, "\n")
stopifnot(n == 1L, gregexpr(old, s, fixed = TRUE)[[1]][1] > 0)
# the mutation: keep the argument, drop the refusal
s <- sub(old, "  if (!ok) invisible(NULL)", s, fixed = TRUE)
writeBin(charToRaw(s), f)
dir.create(MUTLIB, showWarnings = FALSE)
o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(MUTLIB)),
               "--no-multiarch", "--no-docs", "--no-test-load",
               shQuote(MUT)), stdout = TRUE, stderr = TRUE)
cat(tail(o, 2), sep = "\n")
