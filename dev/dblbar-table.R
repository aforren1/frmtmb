# Lane dblbar: the per-case brms comparison table for the findings,
# generated from dev/dblbar-log/compare-{base,lane}.rds and the brms and
# frmtmb probes. Output: dev/dblbar-log/table.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
lg <- function(f) file.path(wt, "dev/dblbar-log", f)
b <- readRDS(lg("compare-base.rds"))
l <- readRDS(lg("compare-lane.rds"))
br <- readRDS(lg("brms-probe.rds"))
fb <- readRDS(lg("frm-probe-base.rds"))
fl <- readRDS(lg("frm-probe-lane.rds"))
yn <- function(x) ifelse(is.na(x), "-", ifelse(x, "yes", "NO"))
short <- function(s) {
  if (is.null(s) || !nzchar(s)) return("")
  substr(gsub("[|]", "/", gsub("\n", " ", s)), 1L, 60L)
}
sink(lg("table.txt"))
cat("| case | spelling | brms sds | frm sds base/lane | Z base/lane |",
    "names base/lane | sd prior rows base/lane | cor base/lane |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(l))) {
  nm <- l$case[i]
  j <- match(nm, b$case)
  sp <- gsub("[|]", "\\\\|", br[[nm]]$case[1])
  if (!is.na(l$Z[i]) || !is.na(b$Z[j])) {
    cat(sprintf("| %s | `%s` | %s | %s/%s | %s/%s | %s/%s | %s/%s | %s/%s |\n",
                nm, sp, l$brms_sds[i], b$frm_sds[j], l$frm_sds[i],
                yn(b$Z[j]), yn(l$Z[i]), yn(b$names[j]), yn(l$names[i]),
                yn(b$prior[j]), yn(l$prior[i]), yn(b$cor[j]), yn(l$cor[i])))
  } else {
    cat(sprintf("| %s | `%s` | brms: %s | base: %s | lane: %s | | | |\n",
                nm, sp, short(br[[nm]]$error %||% "ok"),
                short(fb[[nm]]$error %||% fb[[nm]]$fit_error %||% "ok"),
                short(fl[[nm]]$error %||% fl[[nm]]$fit_error %||% "ok")))
  }
}
ok <- l[!is.na(l$Z), ]
okb <- b[!is.na(b$Z), ]
cat("\nlane: cases compared", nrow(ok), "; Z, names and cor all agree:",
    sum(ok$Z & ok$names & ok$cor), "; sd prior rows agree:", sum(ok$prior),
    "\n")
cat("base: cases compared", nrow(okb), "; Z, names and cor all agree:",
    sum(okb$Z & okb$names & okb$cor), "; sd prior rows agree:",
    sum(okb$prior), "\n")
sink()
cat(readLines(lg("table.txt")), sep = "\n")
