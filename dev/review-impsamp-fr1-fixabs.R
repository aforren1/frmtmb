# A proposed fix as a mutant: a parameter with no column steps by 1e-4 whatever its
# size. A central difference needs no relative step, and an ordinal threshold moves
# the predictor one for one. Source copy outside the worktree.
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
dst <- "C:/Users/adf44/source/r/rev-impsamp-src/fr1-fix_abs"
unlink(dst, recursive = TRUE)
dir.create(file.path(dst, "R"), recursive = TRUE)
dir.create("C:/Users/adf44/source/r/rev-impsamp-lib/fix_abs", recursive = TRUE, showWarnings = FALSE)
for (f in c("DESCRIPTION", "NAMESPACE")) file.copy(file.path(WT, f), dst)
file.copy(list.files(file.path(WT, "R"), full.names = TRUE), file.path(dst, "R"))
for (d in c("inst", "data", "man")) file.copy(file.path(WT, d), dst, recursive = TRUE)
p <- file.path(dst, "R", "importance.R")
s <- readLines(p)
i <- grep("  h <- imp_cross_step * pmax(1, abs(unname(at)))", s, fixed = TRUE)
stopifnot(length(i) == 1L)
s[i] <- "  h <- rep(imp_cross_step, length(at))"
writeLines(s, p)
cat("edited\n")
