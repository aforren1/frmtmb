# Reviewer of lane optrap, item 7: the one-level group with a bound
# (dev/optrap-onelevel.R's design, seeds 101 to 140): what the user is
# told and what the estimates are, per build; and a candidate fix on
# the lane's build, applied in this process only.
#   Rscript dev/optrap-rev-onelevel.R <lib or none> <label> [fix]
#   fix: "none" or "skip1" (no unit start for a block of one level)
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
lab <- a[2]
fix <- if (length(a) > 2) a[3] else "none"
ns <- asNamespace("frmtmb")
if (fix == "skip1") {
  orig <- get("unit_theta_start", ns)
  patched <- function(frame, tpl, units = predictor_units(frame)) {
    # the same, with every block of a one-level grouping left alone
    keep <- vapply(frame[["re_blocks"]], function(bk) {
      nl <- bk[["n_levels"]] %||% bk[["nlev"]] %||%
        tryCatch(nlevels(bk[["group"]]), error = function(e) NA)
      !(is.numeric(nl) && length(nl) == 1L && isTRUE(nl <= 1))
    }, NA)
    fr <- frame
    fr[["re_blocks"]] <- frame[["re_blocks"]][keep]
    orig(fr, tpl, units)
  }
  environment(patched) <- ns
  utils::assignInNamespace("unit_theta_start", patched, "frmtmb")
}
cat(lab, "fix", fix, "lib", dirname(find.package("frmtmb")), "\n")
X <- list()
for (s in 101:140) {
  set.seed(s)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  w <- character()
  f <- withCallingHandlers(
    suppressMessages(frm(bf(y ~ x + (1 | g1)), family = gaussian(),
                         data = d3,
                         prior = set_prior("", class = "b", ub = 0.1))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  e <- f$estimates
  X[[length(X) + 1L]] <- data.frame(
    seed = s, code = f$opt$convergence, obj = f$opt$objective,
    evals = f$opt$evaluations[1], b0 = e$beta[1], b1 = e$beta[2],
    logsig = e$betad[1], logsd = e$theta[1],
    nwarn = length(w),
    warn = substr(gsub("[\t\n]", " ", paste(unique(substr(w, 1, 45)),
                                             collapse = " | ")), 1, 200))
}
X <- do.call(rbind, X)
X$obj <- sprintf("%.10f", X$obj)
write.table(X, sprintf("dev/optrap-rev-log/onelevel-%s.tsv", lab),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(table(code = X$code))
print(X[X$code != 0, c("seed", "code", "obj", "evals", "b0", "b1",
                       "logsig", "logsd", "warn")])
