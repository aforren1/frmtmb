# Review of lane nlstart: is the lane library the worktree's source? Every
# top-level function of R/ parsed and compared with the installed namespace.
source("dev/review-nlstart-common.R")
ns <- asNamespace("frmtmb")
nfun <- 0; ndiff <- 0; nmiss <- 0; diffs <- character(0)
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (!is.call(e) || !as.character(e[[1]]) %in% c("<-", "=")) next
    if (!is.name(e[[2]])) next
    rhs <- e[[3]]
    if (!is.call(rhs) || !identical(rhs[[1]], as.name("function"))) next
    nm <- as.character(e[[2]])
    nfun <- nfun + 1
    if (!exists(nm, envir = ns, inherits = FALSE)) {
      nmiss <- nmiss + 1; next
    }
    g <- get(nm, envir = ns)
    src <- eval(rhs, baseenv())
    if (!identical(deparse(formals(g)), deparse(formals(src))) ||
        !identical(deparse(body(g)), deparse(body(src)))) {
      ndiff <- ndiff + 1; diffs <- c(diffs, nm)
    }
  }
}
cat("functions", nfun, "differ", ndiff, "missing", nmiss, "\n")
if (length(diffs)) cat("differ:", head(diffs, 20), "\n")
for (nm in c("nl_par_units", "nl_body_additive", "nl_sign_patterns",
             "nl_sign_probe", "nl_fit_probed")) {
  cat(nm, "in namespace:", exists(nm, envir = ns, inherits = FALSE), "\n")
}
cat("patterns k=3:\n"); print(do.call(rbind, frmtmb:::nl_sign_patterns(3)))
cat("patterns k=4:\n"); print(do.call(rbind, frmtmb:::nl_sign_patterns(4)))
cat("patterns k=5:\n"); print(do.call(rbind, frmtmb:::nl_sign_patterns(5)))
