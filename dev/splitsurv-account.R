# Every test of the moved files, before (rellib-r11, base tree) and after
# (lane build), matched by test name: where it went and what it asserts.
#   Rscript dev/splitsurv-account.R > dev/splitsurv-out/account.md
rd <- function(dir, pkg, file) {
  f <- file.path(dir, paste0(pkg, "--", sub("[.]R$", "", file), ".tsv"))
  if (!file.exists(f)) return(NULL)
  t <- utils::read.delim(f, stringsAsFactors = FALSE)
  if (!"test" %in% names(t)) return(NULL)
  t$pkg <- pkg
  t
}
B <- "dev/splitsurv-before-out"
A <- "dev/splitsurv-after-out"
moves <- list(
  # before pkg, before file, after pkg, after file
  c("frmtmb.spline", "test-curve.R", "frmtmb", "test-curve.R"),
  c("frmtmb.spline", "test-deriv.R", "frmtmb", "test-curve-deriv.R"),
  c("frmtmb.spline", "test-difference.R", "frmtmb", "test-curve-difference.R"),
  c("frmtmb.spline", "test-gratia.R", "frmtmb", "test-curve-gratia.R"),
  c("frmtmb.spline", "test-illcond.R", "frmtmb", "test-curve-illcond.R"),
  c("frmtmb.spline", "test-new-levels.R", "frmtmb", "test-curve-new-levels.R"),
  c("frmtmb.spline", "test-span.R", "frmtmb", "test-curve-span.R"),
  c("frmtmb.spline", "test-scale.R", "frmtmb", "test-curve-scale.R"),
  c("frmtmb.spline", "test-royston-parmar.R", "frmtmb.survival",
    "test-royston-parmar.R"),
  c("frmtmb.spline", "test-rp-floored.R", "frmtmb.survival",
    "test-rp-floored.R"),
  c("frmtmb.spline", "test-frailty.R", "frmtmb.survival", "test-frailty.R"),
  c("frmtmb.spline", "test-surface.R", "frmtmb.survival",
    "test-rp-surface.R"),
  c("frmtmb.spline", "test-bracket-access.R", "frmtmb.survival",
    "test-bracket-access.R"),
  c("frmtmb.spline", "test-conditions-census.R", "frmtmb.survival",
    "test-conditions-census.R"),
  c("frmtmb.spline", "test-conditions.R", "frmtmb.survival",
    "test-conditions.R"),
  c("frmtmb.spline", "test-message-uniqueness.R", "frmtmb.survival",
    "test-message-uniqueness.R"),
  c("frmtmb", "test-cox-chart.R", "frmtmb.survival", "test-cox-chart.R"),
  c("frmtmb", "test-famgaps.R", "frmtmb", "test-famgaps.R"),
  c("frmtmb", "test-brms-families.R", "frmtmb", "test-brms-families.R"),
  c("frmtmb", "test-hurdle-nb-zoib.R", "frmtmb", "test-hurdle-nb-zoib.R"),
  c("frmtmb", "test-mv-gaps.R", "frmtmb", "test-mv-gaps.R"),
  c("frmtmb", "test-simulate-density.R", "frmtmb", "test-simulate-density.R"),
  c("frmtmb", "test-update-pool.R", "frmtmb", "test-update-pool.R"),
  c("frmtmb.sample", "test-cox-draws.R", "frmtmb.sample", "test-cox-draws.R"),
  c("frmtmb.sample", "test-simulators.R", "frmtmb.sample",
    "test-simulators.R"))
# the other files a moved test may have gone to
extra_after <- list(c("frmtmb", "test-curve.R"),
                    c("frmtmb", "test-simplex-pars.R"),
                    c("frmtmb.survival", "test-cox.R"),
                    c("frmtmb.survival", "test-cox-chart.R"))
cat("| before file | tests | pass before | after file | tests | pass after |\n")
cat("|---|---|---|---|---|---|\n")
for (m in moves) {
  b <- rd(B, m[1], m[2]); a <- rd(A, m[3], m[4])
  cat(sprintf("| %s `%s` | %d | %d | %s `%s` | %d | %d |\n", m[1], m[2],
              if (is.null(b)) 0L else nrow(b),
              if (is.null(b)) 0L else sum(b$passed), m[3], m[4],
              if (is.null(a)) 0L else nrow(a),
              if (is.null(a)) 0L else sum(a$passed)))
}
cat("\nTests of a before file whose name is not in its after file, and where\n")
cat("they are now:\n\n")
pool <- do.call(rbind, c(lapply(moves, function(m) rd(A, m[3], m[4])),
                         lapply(extra_after, function(e) rd(A, e[1], e[2]))))
pool <- unique(pool)
pool$where <- paste(pool$pkg, pool$file)
for (m in moves) {
  b <- rd(B, m[1], m[2]); a <- rd(A, m[3], m[4])
  if (is.null(b)) next
  gone <- setdiff(b$test, if (is.null(a)) character(0) else a$test)
  for (t in gone) {
    hit <- unique(pool$where[pool$test == t])
    pb <- sum(b$passed[b$test == t])
    pa <- sum(pool$passed[pool$test == t])
    cat(sprintf("- %s `%s` \"%s\" (%d passed): %s\n", m[1], m[2], t, pb,
                if (length(hit)) paste0(paste(hit, collapse = "; "),
                                        " (", pa, " passed)") else
                  "NOT FOUND BY NAME"))
  }
}
cat("\nTests in an after file whose name is not in its before file (new or\n")
cat("rewritten):\n\n")
for (m in moves) {
  b <- rd(B, m[1], m[2]); a <- rd(A, m[3], m[4])
  if (is.null(a)) next
  new <- setdiff(a$test, if (is.null(b)) character(0) else b$test)
  for (t in new) {
    cat(sprintf("- %s `%s` \"%s\" (%d passed)\n", m[3], m[4], t,
                sum(a$passed[a$test == t])))
  }
}
for (e in list(c("frmtmb", "test-simplex-pars.R"),
               c("frmtmb", "test-curve-check.R"),
               c("frmtmb", "test-cox-home.R"),
               c("frmtmb", "test-extension-home.R"),
               c("frmtmb.survival", "test-saved-fits.R"),
               c("frmtmb.sample", "test-simplex-family-notes.R"),
               c("frmtmb.survival", "test-cox.R"),
               c("frmtmb.survival", "test-cox-home.R"),
               c("frmtmb.survival", "test-cox-links.R"),
               c("frmtmb.survival", "test-no-reach.R"))) {
  a <- rd(A, e[1], e[2])
  if (is.null(a)) next
  cat(sprintf("- %s `%s`: %d tests, %d passed (new file)\n", e[1], e[2],
              nrow(a), sum(a$passed)))
}
