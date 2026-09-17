# Reviewer, lane wt-priorform: bf() spellings at the edge of the lane's
# change, brms 2.23.0 against the lane.
#   Rscript dev/priorform-rev-bfedge.R brms|ref|lane
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(switch(mode, ref = "C:/Users/adf44/source/r/rellib-r3",
                   lane = "C:/Users/adf44/source/r/priorform-lib", NULL),
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (mode == "brms") suppressMessages(library(brms)) else suppressMessages(library(frmtmb))
show <- function(txt) {
  r <- tryCatch({
    v <- eval(parse(text = txt))
    if (mode == "brms" && inherits(v, "brmsformula"))
      paste("pforms:", paste(names(v$pforms), vapply(v$pforms, function(f) deparse1(f), ""), collapse = "; "),
            "| nl:", isTRUE(attr(v$formula, "nl")))
    else if (inherits(v, "frmtmb_formula"))
      paste("pforms:", paste(names(v$pforms), vapply(v$pforms, function(f) deparse1(f), ""), collapse = "; "),
            "| nl:", v$nl)
    else paste(capture.output(print(v)), collapse = " ")
  }, warning = function(w) paste("WARNING:", conditionMessage(w)),
     error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 110)))
  cat(sprintf("%-55s %s\n", txt, r))
}
for (t in c('bf(y ~ x, phi = sigma ~ z)',
            'bf(y ~ x, sigma = ~ z)',
            'bf(y ~ x, sigma = y ~ z)',
            'bf(bf(y ~ x, sigma ~ z), sigma ~ w)',
            'bf(bf(y ~ x), nl = TRUE)',
            'bf(bf(y ~ a, a ~ 1, nl = TRUE), nl = FALSE)',
            'bf(y ~ x + ~z)',
            'bf(y ~ (~x))',
            'bf(~ ~x)',
            'bf(y ~ x, sigma ~ (~z))',
            'bf(mvbind(y1, y2) ~ a, a ~ 1, nl = TRUE)'))
  show(t)
