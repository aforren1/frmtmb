# Lane illse2, fix round 1: run test-se-reparam.R on the lane library
# with one fix undone in the namespace, so that the tests that pin it
# are seen to fail on the defect.
#   Rscript dev/illse2-fx-mock.R <mode>
# mode:
#   esc17    install 17's escape rule: any optCtrl other than the
#            default turns the sd escape off (review finding 2)
#   schemeE  install 17's difference scheme: Richardson only along a
#            column that moves one parameter (review finding 3)
#   noguard  the asymmetry guard off, se_reparam_asym = Inf (finding 4)
#   allguard the guard refuses every Hessian, se_reparam_asym = 0
#            (finding 4, the absent case of the refusal)
mode <- commandArgs(TRUE)[[1L]]
.libPaths(c("C:/Users/adf44/source/r/wt-illse2-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
put <- function(name, value) {
  unlockBinding(name, ns)
  assign(name, value, envir = ns)
  lockBinding(name, ns)
}
switch(mode,
  esc17 = put("optctrl_stops_earlier", function(oc) {
    !identical(oc, list(iter.max = 1000, eval.max = 1000))
  }),
  schemeE = {
    f <- function(obj, p, M, step = 1e-3, rich_all = FALSE) {
      saved <- obj_state_save(obj)
      on.exit(obj_state_restore(obj, saved), add = TRUE)
      k <- ncol(M)
      H <- matrix(0, k, k)
      cd <- function(j, h) {
        gp <- tryCatch(se_quiet(as.numeric(obj$gr(p + h * M[, j]))),
                       error = function(e) NULL)
        gm <- tryCatch(se_quiet(as.numeric(obj$gr(p - h * M[, j]))),
                       error = function(e) NULL)
        if (length(gp) != length(p) || length(gm) != length(p) ||
              !all(is.finite(gp)) || !all(is.finite(gm))) {
          return(NULL)
        }
        crossprod(M, gp - gm) / (2 * h)
      }
      # install 17: extrapolate only a column that moves one parameter
      rich <- rich_all | colSums(M != 0) == 1L
      for (j in seq_len(k)) {
        c1 <- cd(j, step)
        if (is.null(c1)) return(NULL)
        if (rich[j]) {
          c2 <- cd(j, 2 * step)
          if (is.null(c2)) return(NULL)
          c1 <- (4 * c1 - c2) / 3
        }
        H[, j] <- c1
      }
      list(H = (H + t(H)) / 2, E = abs(H - t(H)) / 2)
    }
    environment(f) <- ns
    put("reparam_hessian", f)
  },
  noguard = put("se_reparam_asym", Inf),
  allguard = put("se_reparam_asym", 0),
  stop("unknown mode ", mode))
cat("mode", mode, "| frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
res <- as.data.frame(test_file(
  "C:/Users/adf44/source/r/frmtmb-wt-illse2/tests/testthat/test-se-reparam.R",
  package = "frmtmb", env = test_env("frmtmb"), reporter = "summary",
  stop_on_failure = FALSE))
print(res[, c("test", "passed", "failed", "error")], row.names = FALSE)
