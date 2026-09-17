.libPaths(c('C:/Users/adf44/source/r/conditions-lib', 'C:/Users/adf44/source/r/rellib-r3', 'C:/Users/adf44/source/r/pinlib', 'C:/Users/adf44/AppData/Local/R/win-library/4.6'))
Sys.setenv(NOT_CRAN = 'true')
suppressPackageStartupMessages(library(frmtmb))
options(warn = 1)
set.seed(1)
res <- tryCatch({
### Name: frmtmb_register_compat
### Title: Contribute to the compatibility matrix from another package
### Aliases: frmtmb_register_compat compat_rule_builder compat_aterm_rules

### ** Examples

# what a contributing package's .onLoad() does
contribute <- function() {
  b <- compat_rule_builder()
  b$r("wiener", "cens()", "refused",
      "The family supplies no lcdf, so there is no CDF to censor with.")
  b$r("wiener", "*", "untested",
      "Not exercised outside this package's own suite.")
  b$r("wiener", "hmm", "untested",
      "Another package's feature, so the rule waits for it.")
  frmtmb_register_compat(features = c("wiener" = "family"),
                         rules = b$rules, expects = "hmm")
}
# the accumulator on its own, which is all a rule set is
b <- compat_rule_builder()
b$r("dec()", "trials()", "refused", "Different response shapes.")
b$rules()



'OK'}, error = function(e) paste('ERROR', paste(class(e), collapse = '/'), conditionMessage(e)))
cat('\nRDEX-RESULT', gsub('\n', ' ', res), '\n')
