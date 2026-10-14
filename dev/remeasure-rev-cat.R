# Reviewer probe: frmtmb's exported categorical family object, for the
# bin-2P categorical log_lik block (tests.log_lik.R:481).
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
fam <- frmtmb::brmsfamily("categorical")
cat("dpars:", fam$dpars, "\n"); print(args(fam$lpdf)); print(body(fam$lpdf))
