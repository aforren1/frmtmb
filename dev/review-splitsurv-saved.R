# Saved RP and cox fits reloaded in a fresh session.
# Usage: Rscript review-splitsurv-saved.R <save|ns|absent>
# save:   lane lib, library(frmtmb.survival), fit and saveRDS
# ns:     lane lib, only library(frmtmb); survival installed, not attached
# absent: rev lib holding only the lane's frmtmb, so survival is not
#         installed (rellib-r11 has frmtmb.spline, not frmtmb.survival)
arm <- commandArgs(trailingOnly = TRUE)[1]
first <- if (arm == "absent") "C:/Users/adf44/source/r/rev-splitsurv-lib" else
  "C:/Users/adf44/source/r/wt-splitsurv-lib"
.libPaths(c(first, "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; survival installed:",
    nzchar(system.file(package = "frmtmb.survival")), "\n")
O <- "C:/Users/adf44/source/r/frmtmb-wt-splitsurv/dev/review-splitsurv-out"
show <- function(lab, expr) {
  r <- tryCatch({
    v <- expr
    if (inherits(v, "frmtmb_fit")) {
      paste0("FIT logLik=", format(as.numeric(logLik(v)), digits = 15))
    } else if (is.data.frame(v)) {
      paste0("OK rows=", nrow(v), " first=", format(v[[".estimate"]][1],
                                                   digits = 10))
    } else paste0("OK ", paste(class(v), collapse = "/"))
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("%-30s %s\n", lab, r))
}
if (arm == "save") {
  suppressMessages(library(frmtmb.survival))
  set.seed(20260905)
  dr <- data.frame(grp = factor(rep(c("A", "B"), each = 300)))
  dr$t <- stats::rweibull(600, shape = 1.6,
                          scale = ifelse(dr$grp == "A", 0.30, 3.0))
  dr$censored <- 0L
  i <- which(dr$grp == "A")[1L]
  dr$t[i] <- 50
  dr$censored[i] <- 1L
  fr <- suppressWarnings(frm(bf(t | cens(censored) ~ grp),
                             family = royston_parmar(df = 3), data = dr))
  bad <- fr
  bad$estimates$betad[] <- bad$estimates$betad - 3
  set.seed(1)
  dd <- data.frame(x = rnorm(150), z = rnorm(150))
  dd$time <- rexp(150, exp(-0.5 + 0.7 * dd$x))
  dd$out <- rbinom(150, 1, 0.3)
  fc <- frm(bf(time | cens(out) ~ x), family = cox(df = 6), data = dd)
  saveRDS(list(fr = fr, bad = bad, fc = fc, dr = dr, dd = dd),
          file.path(O, "saved-fits.rds"))
  cat("saved\n")
}
s <- tryCatch(readRDS(file.path(O, "saved-fits.rds")),
              error = function(e) {
                cat("readRDS ERROR:", conditionMessage(e), "\n"); NULL
              })
cat("survival ns loaded after readRDS:",
    isNamespaceLoaded("frmtmb.survival"), "; attached:",
    "package:frmtmb.survival" %in% search(), "\n")
if (!is.null(s)) {
  gr <- data.frame(t = seq(0.05, 5, length.out = 10),
                   grp = factor("A", levels = levels(s$dr$grp)))
  show("frm_curve(rp good)", frm_curve(s$fr, newdata = gr, dpar = "mu",
                                       simultaneous = FALSE))
  show("frm_curve(rp bad) [gate]", frm_curve(s$bad, newdata = gr,
                                             dpar = "mu",
                                             simultaneous = FALSE))
  show("logLik(rp)", logLik(s$fr))
  show("predict(rp, link)", predict(s$fr, type = "link"))
  show("update(rp, . ~ 1)", suppressWarnings(update(s$fr, . ~ 1)))
  show("update(cox, . ~ . + z)", update(s$fc, . ~ . + z))
  show("summary(cox)", summary(s$fc))
}
