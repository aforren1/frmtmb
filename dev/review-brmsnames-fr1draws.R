## F5: fixef() of frmtmb.sample draws against brms's fixef() of the same
## models (dev/review-brmsnames-out/brms-<n>.rds), cases 8, 13 and 15 of
## dev/review-brmsnames-cases.R, and against the fit's fixef().
.libPaths(c("C:/Users/adf44/source/r/wt-brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(library(frmtmb.sample))
source("dev/review-brmsnames-cases.R")
d <- rv_data()
cases <- rv_cases(frmtmb::bf, frmtmb::mvbf)
for (i in c(8, 13, 15)) {
  ref <- readRDS(sprintf("dev/review-brmsnames-out/brms-%02d.rds", i))
  fit <- q(frm(cases[[i]]$form, data = d))
  ds <- q(frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3))
  fd <- fixef(ds)
  fs <- fixef(ds, summary = FALSE)
  cat("== case", i, cases[[i]]$id, "\n  brms :", ref$fixef,
      "\n  draws:", rownames(fd),
      "\n  draws = brms:", identical(rownames(fd), ref$fixef),
      "; draws = fit:", identical(rownames(fd), rownames(fixef(fit))),
      "; summary = FALSE cols = rows:", identical(colnames(fs),
                                                  rownames(fd)),
      "; nchains attr kept:", !is.null(attr(fs, "nchains")), "\n")
  # each column still holds its own draws after the reordering
  m <- as.matrix(ds)
  ok <- vapply(colnames(fs), function(cn) {
    v <- paste0("b_", cn)
    v %in% colnames(m) && identical(unname(fs[, cn]), unname(m[, v]))
  }, TRUE)
  cat("  columns match their b_ draws:", sum(ok), "of", length(ok), "\n")
}
