# Punch round 1, m1: rebuild the review's four escaping mutants from the
# CURRENT worktree source and install each core into its own library
# under dev/ordinal-p1-mutlib/<id>; dev/ordinal-runtest.R's "mut:<lib>"
# arm then runs the tests against each. Output:
# dev/ordinal-p1-log-mutants.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
out <- file.path(wt, "dev")
muts <- list(
  M2_acat_logit_nodisc = list("R/families.R",
    "      E <- E * disc\n", "      E <- E\n"),
  M2b_grouped_seq_nodisc = list("R/thres.R",
    "    Mj <- disc * (tau[s + pmin(j, nk) - 1L] - eta)",
    "    Mj <- (tau[s + pmin(j, nk) - 1L] - eta)"),
  M2c_cratio_nodisc = list("R/families.R",
    "      M <- M * disc\n", "      M <- M * (if (stopping) disc else 1)\n"),
  M3b_stz_inverse = list("R/thres.R",
    "sum_to_zero = if (ordered) log(diff(tg)) else tg[-k])",
    "sum_to_zero = if (ordered) log(diff(tg)) else tg[-1L])")
)
R <- file.path(R.home("bin"), "R.exe")
Rs <- file.path(R.home("bin"), "Rscript.exe")
keep <- c("DESCRIPTION", "NAMESPACE", "R", "man", "inst", "data",
          "LICENSE", ".Rbuildignore")
sink(file.path(out, "ordinal-p1-log-mutants.txt"), split = TRUE)
for (id in names(muts)) {
  m <- muts[[id]]
  dst <- file.path(out, "ordinal-p1-mut", id)
  lib <- file.path(out, "ordinal-p1-mutlib", id)
  unlink(dst, recursive = TRUE)
  unlink(lib, recursive = TRUE)
  dir.create(file.path(dst, "frmtmb"), recursive = TRUE)
  dir.create(lib, recursive = TRUE)
  for (k in keep) {
    src <- file.path(wt, k)
    if (file.exists(src)) {
      file.copy(src, file.path(dst, "frmtmb"), recursive = TRUE)
    }
  }
  f <- file.path(dst, "frmtmb", m[[1]])
  txt <- paste(readLines(f), collapse = "\n")
  hits <- lengths(regmatches(txt, gregexpr(m[[2]], txt, fixed = TRUE)))
  writeLines(sub(m[[2]], m[[3]], txt, fixed = TRUE), f)
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", lib),
                     "--no-multiarch", "--no-test-load",
                     file.path(dst, "frmtmb")),
                stdout = TRUE, stderr = TRUE)
  cat(id, ": matches", hits, " install:", tail(st, 1), "\n")
}
sink()
