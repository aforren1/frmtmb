# Review of lane superassign: run pkgcheck's own pkgchk_left_assign()
# and its two output functions, sourced verbatim from upstream
# ropensci-review-tools/pkgcheck main (R/checks-left-assign.R, DESCRIPTION
# Version 0.3.4.005), saved as dev/review-superassign-pkgcheck-upstream.R.
# pkgcheck itself is not installed. Its exemptions call pkgstats' ctags
# helpers, which return NULL here when ctags is absent; frmtmb has no
# RefClass and no .onLoad memoise, so no exemption applies either way.
suppressMessages(library(fs))
env <- new.env()
sys.source("C:/Users/adf44/source/r/frmtmb-wt-superassign/dev/review-superassign-pkgcheck-upstream.R",
           envir = env)
R <- "C:/Users/adf44/source/r"
arms <- list(
  base = file.path(R, "frmtmb"),
  lane = file.path(R, "frmtmb-wt-superassign"),
  mutant = file.path(R, "rev-superassign-src/m1"))
pk <- c("", "extensions/frmtmb.ode", "extensions/frmtmb.sample",
        "extensions/frmtmb.survival", "extensions/frmtmb.coupling",
        "extensions/frmtmb.eam", "extensions/frmtmb.latent",
        "extensions/frmtmb.learn", "extensions/frmtmb.spatial")
cat("lib: none (pkgcheck rule sourced from upstream)\n")
for (a in names(arms)) {
  for (p in pk) {
    path <- if (nzchar(p)) file.path(arms[[a]], p) else arms[[a]]
    if (!dir.exists(file.path(path, "R"))) next
    checks <- list(pkg = list(path = path))
    la <- env$pkgchk_left_assign(checks)
    chk <- list(pkg = checks$pkg, checks = list(left_assign = la))
    g <- env$output_pkgchk_global_assign(chk)
    u <- env$output_pkgchk_left_assign(chk)
    # the raw tally of <<- before exemptions, for the count
    fl <- fs::dir_ls(file.path(path, "R"), regexp = "\\.q$|\\.r$|\\.s$",
                     ignore.case = TRUE)
    n <- sum(vapply(fl, function(f) {
      pd <- utils::getParseData(parse(f, keep.source = TRUE))
      sum(pd$token == "LEFT_ASSIGN" & pd$text == "<<-")
    }, 0L))
    nf <- sum(vapply(fl, function(f) {
      pd <- utils::getParseData(parse(f, keep.source = TRUE))
      any(pd$token == "LEFT_ASSIGN" & pd$text == "<<-")
    }, NA))
    cat(sprintf("%-7s %-28s <<- tokens %3d in %2d files; global_assign %s; left_assign %s (<- %d, = %d)\n",
                a, if (nzchar(p)) basename(p) else "frmtmb", n, nf,
                if (g$check_pass) "pass" else "FAIL",
                if (u$check_pass) "pass" else "FAIL",
                la$usage[["<-"]], la$usage[["="]]))
  }
}
