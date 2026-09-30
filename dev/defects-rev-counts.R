# Reviewer of lane defects: one table row per test log, counts summed from
# the RESULT lines; files without a RESULT line are named. Output is
# pasted into dev/reviews/2026-09-29-defects.md.
#   Rscript dev/defects-rev-counts.R > dev/defects-rev-log/counts.md
logs <- c(
  "frmtmb, every file, gates open (lane)" = "core-lane.txt",
  "frmtmb.sample, every file, gates open (lane)" = "sample-lane.txt",
  "frmtmb.latent (lane core)" = "latent-lane.txt",
  "frmtmb.coupling (lane core)" = "ext-frmtmb.coupling.txt",
  "frmtmb.eam (lane core)" = "ext-frmtmb.eam.txt",
  "frmtmb.learn (lane core)" = "ext-frmtmb.learn.txt",
  "frmtmb.ode (lane core)" = "ext-frmtmb.ode.txt",
  "frmtmb.spline (lane core)" = "ext-frmtmb.spline.txt",
  "changed frmtmb files (base, rellib-r3)" = "changed-base.txt",
  "changed frmtmb.sample file (base, rellib-r3)" = "sample-changed-base.txt")
num <- function(r, k) {
  sum(as.integer(sub(paste0(".* ", k, "= ([0-9]+).*"), "\\1", r)))
}
cat("| run | files | pass | fail | error | skip | warn |\n")
cat("|---|---|---|---|---|---|---|\n")
for (k in names(logs)) {
  f <- file.path("dev/defects-rev-log", logs[[k]])
  if (!file.exists(f)) { cat("|", k, "| log missing | | | | | |\n"); next }
  x <- readLines(f)
  r <- grep("^RESULT", x, value = TRUE)
  hdr <- grep("^SUITE", x, value = TRUE)
  nores <- grep("^NO RESULT", x, value = TRUE)
  cat(sprintf("| %s | %s%s | %d | %d | %d | %d | %d |\n", k, hdr,
              if (length(nores)) paste0("; ", paste(nores, collapse = ", "))
              else "", num(r, "pass"), num(r, "fail"), num(r, "err"),
              num(r, "skip"), num(r, "warn")))
}
cat("\nPer file, base and lane, for the files the lane changed or added:\n\n")
cat("| file | base pass/fail/err/skip | lane pass/fail/err/skip |\n")
cat("|---|---|---|\n")
pick <- function(f, name) {
  x <- grep(paste0(" ", name, " "), readLines(f), value = TRUE, fixed = TRUE)
  if (!length(x)) return("not run")
  sprintf("%d/%d/%d/%d", num(x, "pass"), num(x, "fail"), num(x, "err"),
          num(x, "skip"))
}
for (fn in c("test-brms-parity-defects.R", "test-arg-refusal.R",
             "test-compat.R", "test-mv-gaps.R", "test-prior-compat.R")) {
  cat(sprintf("| %s | %s | %s |\n", fn,
              pick("dev/defects-rev-log/changed-base.txt", fn),
              pick("dev/defects-rev-log/core-lane.txt", fn)))
}
cat(sprintf("| frmtmb.sample test-draws-methods.R | %s | %s |\n",
            pick("dev/defects-rev-log/sample-changed-base.txt",
                 "test-draws-methods.R"),
            pick("dev/defects-rev-log/sample-lane.txt",
                 "test-draws-methods.R")))
tier <- c(grep("brms-suite", readLines("dev/defects-rev-log/core-lane.txt"),
               value = TRUE),
          grep("brms-suite", readLines("dev/defects-rev-log/sample-lane.txt"),
               value = TRUE))
cat(sprintf(paste0("\nPorted tier, verdicts asserted (inside the runs above): ",
                   "%d files, pass %d, fail %d, error %d, skip %d, warn %d\n"),
            length(tier), num(tier, "pass"), num(tier, "fail"),
            num(tier, "err"), num(tier, "skip"), num(tier, "warn")))
m <- list.files("dev/defects-rev-log/mutants", full.names = TRUE)
cat("\nMutants (dev/defects-rev-mutant.R):\n\n")
for (f in m) {
  x <- grep("^MUTANT|caught by", readLines(f), value = TRUE)
  cat("-", paste(trimws(x), collapse = "; "), "\n")
}
