# Lane optima: every count dev/optima-findings.md quotes, generated from
# the local logs in dev/optima-log/ (gitignored). Paste its output into
# the findings' generated blocks.
#   Rscript dev/optima-record.R
L <- "dev/optima-log"
rs <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
run <- function(...) {
  out <- system2(rs, c(...), stdout = TRUE, stderr = TRUE)
  cat(out, sep = "\n")
}
rd <- function(p) {
  do.call(rbind, lapply(Sys.glob(file.path(L, paste0(p, "-*.tsv"))),
                        utils::read.delim))
}
cat("### mo(): base (rellib-r6) against the lane, 200 seeds\n")
run("dev/optima-mo-sum.R", file.path(L, "mo-base2"), file.path(L, "mo-lane3"))
cat("\n### mo(): the chart alone (mo_search() disabled) against the lane\n")
run("dev/optima-mo-sum.R", file.path(L, "mo-lanens"), file.path(L, "mo-lane3"))
cat("\n### mo(): plain ls ~ mo(income), 200 seeds\n")
for (a in c("base", "lane3")) {
  X <- utils::read.delim(file.path(L, paste0("mo-main-", a, ".tsv")))
  X$gap <- X$ll_exact - X$ll_fit
  cat(sprintf("  %-5s gap <= 1e-6: %d of %d; max gap %.3g; codes %s; searched %d; evaluations total %d, median %.1f\n",
              a, sum(X$gap <= 1e-6), nrow(X), max(X$gap),
              paste(names(table(X$code)), table(X$code), sep = "=",
                    collapse = " "), sum(X$searched), sum(X$evals),
              stats::median(X$evals)))
}
cat("\n### prototype arms (dev/optima-mo-flip.R)\n")
for (p in c("flip-soft", "flip-sphere", "cone-sphere", "coneskip-sphere",
            "conescreen-sphere")) {
  run("dev/optima-mo-flip-sum.R", file.path(L, p))
}
cat("\n### cs() ordinal mixtures, 20 seeds each\n")
for (a in c("base", "lane6")) {
  X <- rd(paste0("csmix-", a))
  cat("  arm", a, "\n")
  for (dsg in unique(X$design)) {
    tb <- table(X$class[X$design == dsg])
    cat(sprintf("    %-22s %s\n", dsg,
                paste(names(tb), tb, sep = "=", collapse = " ")))
  }
}
A <- rd("csmix-base")
B <- rd("csmix-lane6")
m <- merge(A, B, by = c("design", "seed"))
both <- !is.na(m$ll.x) & !is.na(m$ll.y)
cat("  fits with a logLik on both arms:", sum(both), "; identical:",
    sum(m$ll.x[both] == m$ll.y[both]), "\n")
ctrl <- m$design %in% c("cum", "sratio", "mix_sratio_sratio",
                        "mix_cum_sratio_nocs")
cat("  controls (no cs() on a cumulative component of a mixture):",
    sum(ctrl), "fits; logLik identical", sum(m$ll.x[ctrl] == m$ll.y[ctrl]),
    "\n")
cat("\n### the crossing wall (dev/optima-csmix-wall.R)\n")
cat(utils::tail(readLines(file.path(L, "csmix-wall-lane.txt")), 2), sep = "\n")
cat("\n### ordinal mixture multi-start (dev/optima-ordmix-starts.R)\n")
for (fm in c("cum_cum", "cum_sratio")) {
  X <- rd(paste0("ordstarts-", fm))
  g <- X$gain > 1e-3
  cat(sprintf("  %-10s seeds %d; best of 10 above the default by > 1e-3: %d; of those degenerate: %d; non-degenerate: %d; default degenerate: %d\n",
              fm, nrow(X), sum(g), sum(g & X$degen_best),
              sum(g & !X$degen_best), sum(X$degen0, na.rm = TRUE)))
}
cat("\n### fit_loss2 (dev/optima-loss2.R)\n")
for (f in c("loss2-base.log", "loss2-lane3.log")) {
  x <- readLines(file.path(L, f))
  cat("  ", f, "\n")
  cat(paste0("    ", x[grep("^default|^logLik over|^default below", x)]),
      sep = "\n")
}
