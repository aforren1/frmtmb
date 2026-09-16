# The lane's body guard was seen failing on a mutant that names
# `stancode`. A guard that only catches the one generic its author
# mutated is not a guard, so mutate a DIFFERENT one, in a different
# source file, and check the guard names it.
MUT <- "C:/Users/adf44/source/r/sgrev-mut-loglik"
MUTLIB <- "C:/Users/adf44/source/r/sgrev-mut-loglik-lib"
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen/extensions/frmtmb.sample"
unlink(MUT, recursive = TRUE)
dir.create(MUT, recursive = TRUE)
file.copy(src, MUT, recursive = TRUE)
pkg <- file.path(MUT, "frmtmb.sample")
f <- file.path(pkg, "R", "loo.R")
s <- readChar(f, file.size(f), useBytes = TRUE)
old <- '  UseMethod("log_lik")\n}'
hits <- gregexpr(old, s, fixed = TRUE)[[1]]
cat("the log_lik generic body was found exactly once: ",
    length(hits) == 1L && hits[1] > 0, "\n")
stopifnot(length(hits) == 1L, hits[1] > 0)
s <- sub(old, paste0('  if (inherits(object, "frmtmb_draws")) ',
                     'message("work in the generic")\n',
                     '  UseMethod("log_lik")\n}'), s, fixed = TRUE)
writeBin(charToRaw(s), f)
dir.create(MUTLIB, showWarnings = FALSE)
o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(MUTLIB)),
               "--no-multiarch", "--no-docs", "--no-test-load",
               shQuote(pkg)), stdout = TRUE, stderr = TRUE)
cat(tail(o, 2), sep = "\n")
