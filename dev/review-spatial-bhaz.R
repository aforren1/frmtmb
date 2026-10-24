.libPaths(c("C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/review-spatial-corelib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (fam in c("gaussian()", "'cox'", "brmsfamily('cox')", "cox()"))
  cat(fam, "->", tryCatch({eval(str2lang(sprintf("frm(bf(y | bhaz(df = 5) ~ x), family = %s, data = NULL, dry_run = 'spec')", fam))); "OK"}, error = function(e) conditionMessage(e)), "\n")
