# Does frmtmb_register_compat() refuse a rule naming a feature nobody
# registered yet? (the reason .onLoad() registers mixture_mvn first)
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
r <- tryCatch(frmtmb_register_compat(rules = function() {
  b <- compat_rule_builder()
  b$r("REML", "mixture_mvn", "refused", "probe")
  b$rules()
}), error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
