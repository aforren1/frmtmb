# The body guard seen FAILING: a copy of the fixed frmtmb.sample whose
# `stancode()` generic does work before dispatching, installed into its
# own library, which is then put first for the unmodified test file.
#   Rscript dev/samplegen-mutant.R
MUT <- "C:/Users/adf44/source/r/samplegen-mut-work"
MUTLIB <- "C:/Users/adf44/source/r/samplegen-mut-work-lib"
.libPaths(c("C:/Users/adf44/source/r/samplegen-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
src <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen/extensions/frmtmb.sample"
unlink(MUT, recursive = TRUE)
dir.create(MUT)
file.copy(src, MUT, recursive = TRUE)
pkg <- file.path(MUT, "frmtmb.sample")
f <- file.path(pkg, "R", "methods-draws.R")
s <- readChar(f, file.size(f), useBytes = TRUE)
old <- 'stancode <- function(object, ...) UseMethod("stancode")'
stopifnot(length(gregexpr(old, s, fixed = TRUE)[[1]]) == 1L)
s <- sub(old, paste0('stancode <- function(object, ...) {\n',
                     '  message("work in the generic")\n',
                     '  UseMethod("stancode")\n}'), s, fixed = TRUE)
writeBin(charToRaw(s), f)
dir.create(MUTLIB, showWarnings = FALSE)
o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(MUTLIB)),
               "--no-multiarch", "--no-docs", "--no-test-load",
               shQuote(pkg)), stdout = TRUE, stderr = TRUE)
cat(tail(o, 3), sep = "\n")
