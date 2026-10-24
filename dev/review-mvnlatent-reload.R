# A fit made with library(frmtmb.latent) and an unqualified
# mixture_mvn() call, read back in a session that attaches frmtmb only.
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
show <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    paste("OK", class(v)[1])
  }, error = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("%-45s %s\n", label, substr(r, 1, 300)))
}
s <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-mvnlatent/dev/review-mvnlatent-out/fit3.rds")
fit <- s$fit
dd <- s$dd
cat("latent loaded:", isNamespaceLoaded("frmtmb.latent"), " attached:",
    "package:frmtmb.latent" %in% search(), "\n")
show("predict()", predict(fit))
show("simulate()", simulate(fit, nsim = 1, seed = 2))
show("mixture_probs()", mixture_probs(fit))
cat("latent loaded after methods:", isNamespaceLoaded("frmtmb.latent"), "\n")
show("frm_compat(REML, mixture_mvn)", frm_compat("REML", "mixture_mvn")$status)
show("update(fit)", update(fit))
show("update(fit, newdata = dd)", update(fit, newdata = dd))
