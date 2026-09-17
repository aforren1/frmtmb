## The generated blocks of dev/brmsnames-findings.md, read from the logs
## the measurement scripts wrote. Nothing here is typed from memory.
##
##   Rscript dev/brmsnames-summary.R > dev/brmsnames-log/summary.txt
lg <- function(f) readLines(file.path("dev/brmsnames-log", f), warn = FALSE)
pick <- function(x, pat) grep(pat, x, value = TRUE)

cat("== block 1: the defects, BEFORE (rellib-r3) and AFTER (lane) ==\n")
cat("script dev/brmsnames-probe.R; fit: brms::epilepsy,\n")
cat("count ~ zBase * Trt + (1 | patient), poisson; draws:\n")
cat("dev/brmsnames-draws.R (data seed 9, sampler seed 20260915)\n")
for (arm in c("base", "lane")) {
  x <- lg(sprintf("probe-%s.txt", arm))
  cat("\n-- ", arm, " --\n", sep = "")
  k <- c("^posterior_summary\\(fit\\)", "^names\\(VarCorr", "^class:",
         "^\\[1\\] \"", "^\\[4\\] \"", "^hypothesis\\(fit, 'zBase",
         "^names: ", "^fixef\\(fit, FALSE\\)", "^ranef\\(fit, FALSE\\)",
         "^variables\\(ds\\)", "^bayes_R2\\(ds, NULL, TRUE, TRUE\\)",
         "^posterior_summary\\(ds, ", "^fixef\\(ds, FALSE\\)",
         "^ranef\\(ds, FALSE\\)", "^VarCorr\\(ds, NULL, FALSE\\)")
  for (p in k) {
    for (l in pick(x, p)) cat(strtrim(gsub("\\s+", " ", l), 78), "\n")
  }
}

cat("\n== block 2: names against brms's own parse ==\n")
cat("script dev/brmsnames-naming.R, data seed 41, thirteen models\n")
for (arm in c("base", "lane")) {
  x <- lg(sprintf("naming-%s.txt", arm))
  cat(pick(x, "^== total"), "\n")
  lacks <- trimws(sub("^names brms lacks: *", "",
                      pick(x, "^names brms lacks")))
  cat(sprintf("  %s: models with a name brms lacks: %d of %d\n", arm,
              sum(nzchar(lacks)), length(lacks)))
}
cat("models:\n")
cat(paste0("  ", sub("^== (.*) ==$", "\\1",
                     pick(lg("naming-lane.txt"), "^== [^t]"))), sep = "\n")

cat("\n== block 3: output match against brms's installed methods ==\n")
cat("script dev/brmsnames-match.R; 45 calls on each of two draws sets\n")
for (arm in c("base", "lane")) {
  tab <- utils::read.delim(sprintf("dev/brmsnames-log/match-%s.tsv", arm),
                           stringsAsFactors = FALSE, quote = "")
  ctl <- startsWith(tab$call, "CONTROL")
  cat(sprintf("%-5s comparisons %d:", arm, sum(!ctl)))
  tt <- table(tab$result[!ctl])
  cat(paste0(" ", names(tt), " ", as.integer(tt)), "\n")
  cc <- tab[ctl, ]
  cat(sprintf("      control %-50s %s\n", substr(cc$call, 9, 58),
              cc$result), sep = "")
  if (arm == "lane") {
    nid <- tab[!ctl & tab$result != "identical", c("model", "call",
                                                   "result")]
    cat("      not identical():",
        paste(nid$model, nid$call, nid$result), sep = "\n        ")
    cat("\n")
  }
}

cat("\n== block 4: positional formals audit of 68 draws methods ==\n")
cat("script dev/brmsnames-beyond.R (dev/brmsmatch-beyond.R's criterion)\n")
for (arm in c("base", "lane")) {
  x <- lg(sprintf("beyond-%s.txt", arm))
  cat(arm, ": ", gsub("\\s+", " ", pick(x, "^diverge:")), "\n", sep = "")
  cat(arm, ": ", gsub("\\s+", " ", pick(x, "^scored 'agree'")), "\n",
      sep = "")
}
cat("lane, extra positional argument on the methods still shorter\n")
cat("than brms's (dev/brmsnames-dotcheck.R):\n")
x <- lg("dotcheck-lane.txt")
x <- x[x != "DONE"]
cat("  refused:", sum(grepl("refused:", x)), " answered:",
    sum(grepl("ANSWERED", x)), "\n")
cat(paste0("  ", pick(x, "ANSWERED")), sep = "\n")

cat("\n== block 5: blast radius, lines reading a changed surface ==\n")
cat("script dev/brmsnames-blast.R, run on the unchanged tree\n")
x <- lg("blast-base.txt")
i <- grep("^== lines per surface, by kind", x)
j <- grep("^== lines per surface, by package", x)
cat(x[(i + 1):(j - 1)], sep = "\n")

cat("\n== block 6: test files, one per R process ==\n")
cat("runner dev/brmsnames-runtest.R, sums failures AND errors\n")
num_after <- function(x, k) {
  as.integer(vapply(regmatches(x, regexec(paste0(k, " +([0-9]+)"), x)),
                    function(m) m[2L], ""))
}
tally <- function(lines) {
  x <- grep("^BLOCKS", lines, value = TRUE)
  sprintf("%3d files %4d blocks PASS %5d FAIL %d ERROR %d SKIP %d",
          length(x), sum(num_after(x, "BLOCKS")), sum(num_after(x, "PASS")),
          sum(num_after(x, "FAIL")), sum(num_after(x, "ERROR")),
          sum(num_after(x, "SKIP")))
}
cat("core, full suite:    ", tally(lg("full-core-summary.txt")), "\n")
ext <- lg("full-ext-summary.txt")
for (p in c("coupling", "eam", "latent", "learn", "ode", "sample",
            "spline")) {
  cat(sprintf("%-20s ", paste0(p, ", full suite:")),
      tally(grep(paste0("frmtmb[.]", p, "/"), ext, value = TRUE)), "\n")
}
one <- function(f) {
  if (!file.exists(f)) return("(not run)")
  x <- readLines(f, warn = FALSE)
  g <- function(k) sum(num_after(grep(paste0("^", k, " "), x,
                                      value = TRUE), k))
  sprintf("PASS %d FAIL %d ERROR %d", g("PASS"), g("FAIL"), g("ERROR"))
}
cat("\nthe new and changed pinning tests, base build then lane build\n")
cat("(the lane column reads the full-suite logs above):\n")
d <- "dev/brmsnames-log"
rows <- list(
  c("test-brms-names.R", "newtests-base", "full-core",
    "tests_testthat_test-brms-names.R.log"),
  c("test-brms-output.R", "newtests-base", "full-ext",
    "extensions_frmtmb.sample_tests_testthat_test-brms-output.R.log"),
  c("test-brms-pins.R", "newtests-base", "full-ext",
    "extensions_frmtmb.sample_tests_testthat_test-brms-pins.R.log"),
  c("test-counterfactual.R", "learn-base", "full-ext",
    "extensions_frmtmb.learn_tests_testthat_test-counterfactual.R.log"))
for (r in rows) {
  cat(sprintf("  %s\n    base: %s\n    lane: %s\n", r[1L],
              one(file.path(d, r[2L], r[4L])),
              one(file.path(d, r[3L], r[4L]))))
}
cat("\nwith an EMPTY Stan cache (dev/brmsnames-runtest-nocache.R):\n")
for (f in c("nocache-brms-output.log", "nocache-brms-pins.log",
            "nocache-brms-methods.log")) {
  cat(sprintf("  %-28s %s\n", f, one(file.path(d, f))))
}

cat("\n== block 7: punch rounds 1 and 2, the pins under mutation ==\n")
cat("script dev/brmsnames-mutants.R <mutant> on test-brms-pins.R,\n")
cat("test-brms-output.R, test-brms-names.R and frmtmb.latent's\n")
cat("test-hmm.R (lane build, in memory); from mixnat on, round 2\n")
for (m in c("none", "rlev", "rcoef", "bcoef", "corord", "natsig",
            "residse", "llpw", "norename", "charmap", "stanname", "lvldot",
            "lvljoin", "nodupref", "nosuffix", "noredup", "natname",
            "rrfill", "lapref", "mvresid", "mvr2", "barenp", "mixnat",
            "hmmnat", "rdup", "lvlmerge", "respdup", "hypdup", "sdscolon",
            "bspname")) {
  x <- lg(file.path("mutants", paste0(m, ".txt")))
  r <- pick(x, "^test-")
  cat(sprintf("  %-9s %s\n", m,
              paste(gsub("\\s+", " ", sub("^test-(brms-)?", "", r)),
                    collapse = " | ")))
}

cat("\n== block 8: punch round 1, BLOCKER, the colon pin as a script ==\n")
cat("script dev/brmsnames-pin-colon.R, data seed 8\n")
for (arm in c("base", "lane")) {
  x <- lg(sprintf("pin-colon-%s.txt", arm))
  cat(paste0("  ", arm, ": ", pick(x, "^hypothesis|^pin:|^b_x =")),
      sep = "\n")
}

cat("\n== block 9: punch round 1, collisions measured on brms ==\n")
cat("brms: dev/brmsnames-rev-log/collide-brms.txt (reviewer's script);\n")
cat("lane: dev/brmsnames-rev-collide.R lane, data seed 12\n")
x <- lg("collide-lane.txt")
cat(paste0("  ", pick(x, paste0("^==|^variables:|^frm ERROR|",
                                "^draws labels dup|^hypothesis\\("))),
    sep = "\n")

cat("\n== block 10: round 1 verified on the current build ==\n")
cat("MAJOR 2, which coefficient takes __1, against the reviewer's C6\n")
cat("brmsfit (dev/brmsnames-verify-suffix.R):\n")
cat(paste0("  ", lg("verify-suffix.txt")), sep = "\n")
cat("MAJOR 4, laplace draws (dev/brmsnames-probe-laplace.R, data seed 5,\n")
cat("sampler seed 8):\n")
cat(paste0("  ", strtrim(lg("probe-laplace.txt"), 76)), sep = "\n")

cat("\n== block 11: punch round 2 ==\n")
cat("BLOCKER, mixture weights: dev/brmsnames-rev2-mixture.R lane\n")
cat("(reviewer's script, data seed 64). brms: theta1 0.8759, theta2\n")
cat("0.1241, hypothesis(theta1 = 0.5) 0.3759\n")
x <- lg("rev2-mixture-lane.txt")
cat(paste0("  ", strtrim(x[grep("^variables|^1 |^theta1 |ERROR|theta2 *$",
                                x)], 76)), sep = "\n")
i <- grep("^== draws column means", x)
cat(paste0("  ", x[i + 1:4]), sep = "\n")
cat("MAJOR and MINORs, collisions: dev/brmsnames-rev2-collide.R lane\n")
cat("(data seed 52); C1 is refused, so its section prints nothing\n")
x <- lg("rev2-collide-lane.txt")
cat(paste0("  ", strtrim(pick(x, "^==|frm |draws labels|as_draws_df|sds"),
                         76)), sep = "\n")
cat("natural-flag audit (dev/brmsnames-natural-audit.R): every non-primary\n")
cat("dpar of every family in core, eam, latent and learn can take the\n")
cat("flag; the rows not read elementwise, and the count of the rest:\n")
x <- lg("natural-audit.txt")
x <- x[seq_len(grep("^joint rows", x) - 1L)]
cat(paste0("  ", pick(x, "softmax")), sep = "\n")
cat(sprintf("  elementwise rows, natural through linkinv: %d\n",
            length(pick(x, "natural, linkinv"))))
cat("DONE\n")
