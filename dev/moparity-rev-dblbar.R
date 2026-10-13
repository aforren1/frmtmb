# Reviewer: || with a factor, with and without mo(), frmtmb vs brms names.
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(4)
n <- 400
dd <- data.frame(x = sample(0:3, n, TRUE), z = rnorm(n),
                 f2 = factor(sample(c("a", "b"), n, TRUE)),
                 g = factor(sample(paste0("g", 1:25), n, TRUE)))
dd$y <- 1 + 0.5 * dd$x + rnorm(25)[dd$g] + rnorm(n)
for (fo in list(y ~ f2 + (f2 || g), y ~ mo(x) * f2 + (mo(x) * f2 || g),
                y ~ mo(x) * f2 + (mo(x) * f2 | g))) {
  r <- tryCatch(suppressWarnings(frm(bf(fo), family = gaussian(), data = dd)),
                error = function(e) e)
  cat("frmtmb", deparse1(fo), ":", if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1, 100)) else
    paste(vapply(r$frame$re_blocks, function(b) paste(unlist(lapply(b$components, `[[`, "cnms")), collapse = ","), ""), collapse = " | "), "\n")
  sc <- brms::stancode(fo, data = dd)
  cat("brms  ", deparse1(fo), ":", paste(trimws(grep("mu[[]n[]] [+]=", strsplit(sc, "\n")[[1]], value = TRUE)), collapse = " ;; "), "\n")
}
