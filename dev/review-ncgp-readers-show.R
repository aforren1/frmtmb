source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
a <- readRDS(file.path(REV_OUT, "readers-lane.rds"))
b <- readRDS(file.path(REV_OUT, "readers-base.rds"))
cat("== ll, conv, theta\n")
for (nm in names(a)) cat(sprintf("%-8s ll %.5f %.5f conv %d %d theta %s | %s\n", nm, a[[nm]]$ll, b[[nm]]$ll, a[[nm]]$conv, b[[nm]]$conv,
    paste(round(a[[nm]]$theta, 3), collapse = ","), paste(round(b[[nm]]$theta, 3), collapse = ",")))
cat("== refusals\n")
ref <- list()
for (nm in names(a)) for (k in names(a[[nm]])) {
  x <- a[[nm]][[k]]
  if (is.character(x) && length(x) == 1 && grepl("^ERR", x)) ref[[length(ref) + 1]] <- data.frame(k = k, m = nm, x = substr(x, 1, 100))
}
ref <- do.call(rbind, ref)
print(ref[!duplicated(ref[, c("k", "x")]), ], row.names = FALSE)
show <- function(nm, k) { cat("--", nm, k, "\nlane\n"); print(a[[nm]][[k]]); cat("base\n"); print(b[[nm]][[k]]) }
show("gp", "emtr"); show("gp", "deriv"); show("gp", "deriv2"); show("byw", "deriv"); show("re", "deriv2")
show("two", "theta"); show("sigma", "theta"); show("grF", "mfx"); show("binom", "mfx")
