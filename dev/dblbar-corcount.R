# Lane dblbar: how many probe cases have correlation lines in brms's
# Stan code. Reads dev/dblbar-log/brms-probe.rds.
br <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-dblbar/dev/dblbar-log/brms-probe.rds")
ok <- Filter(function(b) is.null(b$error), br)
n_cor <- vapply(ok, function(b) length(b$corlines), 1L)
cat("cases with brms Stan code:", length(ok), "\n")
cat("with no cor_/L_/Cor_ line:", sum(n_cor == 0L), "\n")
cat("with one:", paste(names(n_cor)[n_cor > 0L], collapse = ", "), "\n")
