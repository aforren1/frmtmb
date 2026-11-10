# Reviewer scratch: frm()'s message on each of the 40 pairs, plus
# extension aterms on an ordinal family, plus non-ordinal absent cases.
# Seed 20261110.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), find.package("frmtmb"),
    "\n")
set.seed(20261110)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), nt = 4L, cc = 0L, s = 1,
                k = 1L, rr = 0.5, tt = 2, w = 1,
                two = rep(c("upper", "lower"), length.out = n))
d$y <- as.integer(cut(0.8 * d$x + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf)))
spell <- c(trials = "trials(nt)", cens = "cens(cc)", trunc = "trunc(lb = 0)",
           se = "se(s)", mi = "mi()", vint = "vint(k)", vreal = "vreal(rr)",
           rate = "rate(tt)")
mk <- function(at, st) {
  lhs <- paste0("y | ", at)
  switch(st,
    thres = list(f = paste0("y | thres(3) + ", at, " ~ x"), fam = cumulative()),
    disc = list(f = c(paste(lhs, "~ x"), "disc ~ 0 + z"), fam = cumulative()),
    equidistant = list(f = paste(lhs, "~ x"),
                       fam = cumulative(threshold = "equidistant")),
    sum_to_zero = list(f = paste(lhs, "~ x"),
                       fam = cumulative(threshold = "sum_to_zero")),
    cs = list(f = paste(lhs, "~ x + cs(z)"), fam = sratio()),
    plain = list(f = paste(lhs, "~ x"), fam = cumulative()))
}
run <- function(at, st, fam = NULL) {
  m <- mk(at, st)
  dd <- d
  if (grepl("^mi", at)) dd$y[c(3, 9)] <- NA
  form <- do.call(bf, lapply(m$f, stats::as.formula))
  r <- tryCatch({
    suppressWarnings(suppressMessages(frm(form, family = m$fam, data = dd)))
    "FITS"
  }, error = function(e) paste("ERR:", conditionMessage(e)))
  substr(gsub("\n", " ", r), 1, 200)
}
for (at in names(spell)) for (st in c("plain", "thres", "disc", "equidistant",
                                       "sum_to_zero", "cs")) {
  cat(sprintf("%-6s %-12s %s\n", at, st, run(spell[[at]], st)))
}
for (p in c("frmtmb.eam", "frmtmb.survival")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("\n-- extension aterms on an ordinal family\n")
cat("dec x thres :", run("dec(two)", "thres"), "\n")
cat("dec x cs    :", run("dec(two)", "cs"), "\n")
cat("bhaz x thres:", run("bhaz()", "thres"), "\n")
cat("bhaz x cs   :", run("bhaz()", "cs"), "\n")
cat("frm_compat dec() x thres():", frm_compat("dec()", "thres()")$status,
    "\n")
cat("frm_compat bhaz() x cs_pred():", frm_compat("bhaz()", "cs_pred()")$status,
    "\n")
