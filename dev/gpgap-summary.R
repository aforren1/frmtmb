# The findings' number blocks, generated from the logs in dev/gpgap-out/
# and dev/gpgap-equiv-out/, so no count in dev/gpgap-findings.md is typed.
#   Rscript dev/gpgap-summary.R
L <- function(f) readLines(f, warn = FALSE)
cat("### perm (dev/gpgap-perm.R, base build, before the fix; p094 printed to the terminal: A -94.434083, A1 -0.028199)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/perm-[0-9].txt"))) {
  x <- L(f)
  cat(x[grepl("GPPERM|  A|sorted|off-diagonal", x)], sep = "\n")
}
cat("```\n\n### nugget line (dev/gpgap-nugget.R)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/nug-*.txt"))) {
  x <- L(f); cat(x[grepl("GPNUG|constant|gradient", x)], sep = "\n")
}
cat("```\n\n### extended precision (dev/gpgap-mpfr.R)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/mpfr-*.txt"))) {
  x <- L(f); i <- grep("GPMPFR", x); cat(x[i:length(x)], sep = "\n")
}
cat("```\n\n### the review's numbers (dev/gpgap-review83.R)\n\n```\n")
x <- L("dev/gpgap-out/review83.txt"); cat(x[grepl("GP83", x)], sep = "\n")
cat("```\n\n### identities (dev/gpgap-identity.R, seed 20261102)\n\n```\n")
for (i in 1:15) {
  x <- L(sprintf("dev/gpgap-out/id-%d.txt", i)); cat(x[grepl("GPID", x)], sep = "\n")
}
cat("```\n\n### nugget variants (dev/gpgap-nugvar.R)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/var-*.txt"))) {
  x <- L(f)
  cat(x[grepl("GPVAR|^  rel|^  abs|closed form|brms's model", x)], sep = "\n")
}
cat("```\n\n### low noise (dev/gpgap-lownoise.R, seed 20261102)\n\n```\n")
x <- L("dev/gpgap-out/lownoise.txt"); cat(x[grepl("GPLOW|ncp|k = 30", x)], sep = "\n")
cat("```\n\n### equiv pairs against the baseline\n\n")
base <- "C:/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out"
cat("| pair | baseline verdict (level) | L1 const | now (level) | L1 const |",
    "at the 2 points | max grad | brms's own line | s |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (p in c("p022", "p023", "p094", "p109", "p136")) {
  a <- utils::read.delim(file.path(base, paste0(p, ".tsv")))
  b <- utils::read.delim(file.path("dev/gpgap-equiv-out", paste0(p, ".tsv")))
  cat(sprintf("| %s | %s (%s) | %.4g | %s (%s) | %.3g | %s | %.2g | %s | %.0f |\n",
              p, a$verdict, a$level, a$l1_const, b$verdict, b$level, b$l1_const,
              b$l1_extra, b$l1_grad, sub("^brms.s own line reads (\\S+) .*", "\\1", b$l1_by_decision),
              b$seconds))
}
cat("\n### the warning's calibration (dev/gpgap-calib.R, dev/gpgap-calib-sum.R)\n\n```\n")
x <- L("dev/gpgap-out/calib-sum.txt")
cat(x[grepl("^rows|^old|^r_|^scan|^fix", x)], sep = "\n")
cat("the rule is r_q90\n")
rows <- x[grepl("TRUE|FALSE", x) & !grepl("^old|^r_", x)]
tk <- strsplit(trimws(rows), " +")
r <- vapply(tk, function(t) t[length(t)], "")
mv <- vapply(tk, function(t) t[length(t) - 2L], "")
bg <- vapply(tk, function(t) t[length(t) - 1L], "")
cat("MISS:\n"); cat(rows[(mv == "TRUE" | bg == "TRUE") & r == "FALSE"], sep = "\n")
cat("FALSE ALARM:\n"); cat(rows[mv == "FALSE" & bg == "FALSE" & r == "TRUE"], sep = "\n")
fl <- grepl("^ *(floor|scale|xs10)", rows)
cat(sprintf("one row per position (floor, scale, xs10): %d of %d fire, %d of those below 1 SE\n", sum(fl & r == "TRUE"), sum(fl), sum(fl & r == "TRUE" & mv == "FALSE" & bg == "FALSE")))
cat("```\n\n### binomial at 1e5 trials (dev/gpgap-binom1e5.R)\n\n```\n")
x <- L("dev/gpgap-out/binom1e5.txt"); cat(x[grepl("ll |residual|nugget variance", x)], sep = "\n")
cat("```\n\n### mutants of the warning (dev/gpgap-mut.R)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/mut/*.txt"))) { x <- L(f); cat(x[grepl("^MUT", x)], sep = "\n") }
cat("```\n\n### the sampling route (dev/gpgap-sampler.R, dev/gpgap-sampler-l1.R, dev/gpgap-marg15.R)\n\n```\n")
for (f in sort(Sys.glob("dev/gpgap-out/samp/l1-*.txt"))) { x <- L(f); cat(x[grepl("^GPSAMPL1", x)], sep = "\n") }
x <- L("dev/gpgap-out/samp/marg15.txt"); cat(x[grepl("^MARG", x)], sep = "\n")
cat("-- round 2, the installed build, farthest-first factor\n")
for (f in sort(Sys.glob("dev/gpgap-out/samp/r2-l1-*.txt"))) { x <- L(f); cat(x[grepl("^GPSAMPL1", x)], sep = "\n") }
for (f in c("r2-low-frm-15", "r2-low-frm-40", "r2-low-fit-15", "r2-low-fitlap-15", "r2-low-lap-15", "r2-low-lap-40", "r2-units-frm-60", "r2-units-lap-60")) {
  x <- L(file.path("dev/gpgap-out/samp", paste0(f, ".txt"))); cat(x[grepl("^GPSAMP|^frm[(][)] sigma", x)], sep = "\n")
}
x <- L("dev/gpgap-out/r2/ordercheck.txt"); cat(x[grepl("^GPORD", x)], sep = "\n")
cat("lapcost, now, then the centered route at 1e-6 (dev/review-gpgap-lapcost.R, one run):\n")
for (f in c("lapcost-now", "lapcost-base")) { x <- L(file.path("dev/gpgap-out/r2", paste0(f, ".txt"))); cat(x[grepl("^LCOST", x)], sep = "\n") }
for (f in sort(Sys.glob("dev/gpgap-out/mut2/*.txt"))) { x <- L(f); cat(x[grepl("^MUT2", x)], sep = "\n") }
cat("-- before: the factor over sorted positions (fix round 1)\n")
for (f in c("low-lap-15", "low-baselap-15", "low-lap-40", "low-baselap-40", "units-lap-60", "low-frm-15-diag", "low-frm-15", "low-brms-15", "units-frm-60", "units-brms-60", "try-frm-10", "try-brms-10", "try-frm-15-dense", "try-brms-15-dense")) {
  x <- L(file.path("dev/gpgap-out/samp", paste0(f, ".txt"))); cat(x[grepl("^GPSAMP", x)], sep = "\n")
}
cat("```\n\n### k and the rest\n\n```\n")
x <- L("dev/gpgap-out/kstart.txt"); cat(x[grepl("^KS", x)], sep = "\n")
x <- L("dev/gpgap-out/msgcost.txt"); cat(x[grepl("^GPCOST", x)], sep = "\n")
cat("```\n")
cat("\n### test files\n\n")
for (d in sort(Sys.glob("dev/gpgap-out/*-gated"))) {
  fs <- Sys.glob(file.path(d, "*.txt"))
  pk <- sub("--.*", "", basename(fs))
  for (k in unique(pk)) {
    res <- unlist(lapply(fs[pk == k], function(f) {
      x <- grep("^RESULT", readLines(f, warn = FALSE), value = TRUE)
      if (length(x)) x[length(x)] else paste("NO RESULT", basename(f))
    }))
    ok <- grepl("^RESULT", res)
    num <- function(w) {
      sum(as.integer(sub(paste0(".*", w, "=([0-9]+).*"), "\\1", res[ok])))
    }
    libs <- unique(unlist(lapply(fs[pk == k], function(f) {
      x <- grep("^(lib|core lib):", readLines(f, warn = FALSE), value = TRUE)
      sub(" version.*", "", x)
    })))
    cat(sprintf("- %s (%s): %d files, pass %d, fail %d, error %d, skip %d, warn %d. %s\n",
                k, basename(d), length(res), num("pass"), num("fail"),
                num("error"), num("skip"), num("warn"),
                paste(libs, collapse = "; ")))
    bad <- res[!grepl("fail=0 error=0", res) | grepl("warn=[1-9]", res)]
    if (length(bad)) cat(paste0("  - ", bad, "\n"), sep = "")
  }
}
cat("\n### OpenBLAS, LAPACK routed (dev/ciharden-openblas.sh); each log names its BLAS\n\n")
for (v in c("26", "32")) for (f in sort(Sys.glob(sprintf("dev/gpgap-out/ob%s/*.txt", v)))) {
  x <- L(f)
  cat(sprintf("- 0.3.%s: `%s`; `%s`\n", v, x[grepl("^RESULT", x)][1], substr(x[grepl("^blas", x)][1], 1, 90)))
}
if (file.exists("dev/gpgap-out/gpgapsess.log")) {
  x <- L("dev/gpgap-out/gpgapsess.log")
  cat("\n### one session (ungated, nine packages, three orders)\n\n```\n")
  cat(x[grepl("^SESSION", x)], sep = "\n")
  cat("```\n")
  cat(sprintf("\nLEAK lines %d, untagged %d.\n", sum(grepl("LEAK", x)), sum(grepl("LEAK", x) & !grepl("[[]with load", x))))
}
for (f in c("frmtmb", "frmtmb.spatial", "frmtmb.sample")) {
  g <- file.path("dev/gpgap-out", paste0("check-", f, ".txt"))
  if (file.exists(g)) cat(sprintf("\n`R CMD check --as-cran` %s: `%s`\n", f, L(g)[1]))
}
