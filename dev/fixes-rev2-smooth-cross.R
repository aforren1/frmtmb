# Reviewer of lane fixes, re-check: for the fits where one build ends
# below the other, start the worse build at the better build's smoothing
# SDs and dispersion (the random part spans the same covariance in both
# bases, so theta and betad are comparable) and see whether it reaches
# the better logLik: then the model is the same and the worse build
# stopped early.
#   Rscript dev/fixes-rev2-smooth-cross.R <lib of the worse build> <other tag>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
src <- readLines("dev/fixes-rev2-smooth.R")
eval(parse(text = src[grep("^gs <- function", src):
                        (grep("^rows <- list", src) - 1L)]))
other <- {
  f <- list.files("dev/fixes-rev2-log", paste0("^smooth-", a[2], "-.*rds$"),
                  full.names = TRUE)
  x <- unlist(lapply(f, readRDS), recursive = FALSE)
  names(x) <- vapply(x, function(r) paste(r$case, r$seed), "")
  x
}
for (k in strsplit(a[3], ",")[[1]]) {
  cs <- cases[[sub(" .*", "", k)]]
  seed <- as.integer(sub(".* ", "", k))
  d <- if (cs[[1]] == 2) gs(2, seed, n = 300) else gs(cs[[1]], seed)
  f0 <- suppressWarnings(frm(cs[[2]], data = d, family = cs[[3]]))
  o <- other[[k]]
  tpl <- par_template(cs[[2]], data = d, family = cs[[3]])
  st <- list(theta = unname(o$theta))
  ob <- f0$opt$par[names(f0$opt$par) == "betad"]
  fit <- tryCatch(suppressWarnings(frm(cs[[2]], data = d, family = cs[[3]],
                                       start = st)),
                  error = function(e) conditionMessage(e))
  cat(sprintf("%-7s default start ll %.6f | started at %s's theta: %s | %s ll %.6f\n",
              k, as.numeric(logLik(f0)), a[2],
              if (is.character(fit)) fit else sprintf("ll %.6f conv %d",
                                                      as.numeric(logLik(fit)),
                                                      fit$opt$convergence),
              a[2], o$ll))
}
