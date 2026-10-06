# Lane fixes, punch round: dev/fixes-sx-conv4.R under frmtmb_control(autoscale = TRUE).
# fixef() standard errors are not finite on the lane build: the
# smallest smooth log-sd on each build, and what the lane build says.
#   Rscript dev/fixes-sx-conv4.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
forms <- list(y ~ s(x1) + s(x2), y ~ s(x1, bs = "cr", k = 6),
              y ~ s(x1, by = g) + g, y ~ s(x1, by = z), y ~ s(x1, x2),
              y ~ t2(x1, x2), y ~ s(x0) + s(x1) + s(x2) + s(x3))
cases <- rbind(c(1, 3), c(2, 1), c(2, 7), c(3, 1), c(3, 3), c(3, 7),
               c(4, 7), c(6, 3), c(11, 7), c(15, 1), c(17, 1), c(17, 3),
               c(17, 7), c(18, 7), c(19, 1), c(20, 1), c(20, 3), c(20, 7),
               c(13, 7))
for (i in seq_len(nrow(cases))) {
  set.seed(cases[i, 1])
  d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                     verbose = FALSE))
  d$z <- runif(200)
  d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
  fit <- suppressWarnings(frm(bf(forms[[cases[i, 2]]]), data = d, control = frmtmb_control(autoscale = TRUE)))
  th <- fit$estimates$theta
  w <- character()
  se <- withCallingHandlers(fixef(fit)[, "Est.Error"],
                            warning = function(x) {
                              w <<- c(w, conditionMessage(x))
                              invokeRestart("muffleWarning")
                            })
  cat(sprintf("seed %2d form %d: min theta %7.2f  finite SE %s  %s\n",
              cases[i, 1], cases[i, 2], min(th), all(is.finite(se)),
              if (length(w)) substr(w[1], 1, 110) else "(no warning)"))
}
