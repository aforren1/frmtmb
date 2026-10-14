# Reviewer round 2: spot checks of review items 2, 3 and the small ones.
# Output: dev/famsimplex-rev2-spot.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  r <- tryCatch(suppressMessages(suppressWarnings({
    v <- expr
    paste("OK:", paste(utils::capture.output(print(v)), collapse = " | "))
  })), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("\n---", lab, "\n", r, "\n")
}
dd <- fs_data_dirichlet(seed = 11)
dl <- fs_data_logistic_normal(seed = 13)
d2 <- fs_data_dirichlet2(seed = 15)
fits <- list(
  dirichlet = frm(bf(Y ~ x), family = dirichlet(), data = dd),
  logistic_normal = frm(bf(Y ~ x), family = logistic_normal(), data = dl),
  dirichlet2 = frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = d2))
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n==========", nm, "\n")
  show("dharma_residuals", dharma_residuals(f, nsim = 20, seed = 1))
  for (ty in c("response", "pearson", "deviance", "osa")) {
    show(paste("residuals", ty), residuals(f, type = ty))
  }
  show("residuals newdata", residuals(f, newdata = f$frame$data_frame[1:3, ]))
  show("CE method predict", conditional_effects(f, "x",
                                                method = "posterior_predict"))
  show("CE method predict, dpar suggested", conditional_effects(f, "x",
       dpar = sub("^mu", "mu", f$spec$responses[[1]]$primary_dpars[1])))
}
show("dirichlet(refcat = 2), unnamed", {
  d <- dd
  d$Y <- unname(dd$Y)
  rownames(fixef(frm(bf(Y ~ x), family = dirichlet(refcat = 2), data = d)))
})
show("dirichlet(refcat = 2), named a b c", frm(bf(Y ~ x),
     family = dirichlet(refcat = 2), data = dd))
show("logistic_normal(refcat = 2), unnamed", {
  d <- dl
  d$Y <- unname(dl$Y)
  rownames(summary(frm(bf(Y ~ x), family = logistic_normal(refcat = 2),
                       data = d))$spec_pars)
})
show("dirichlet(refcat = 2.5)", dirichlet(refcat = 2.5)$simplex)
show("dirichlet(refcat = TRUE)", dirichlet(refcat = TRUE))
show("brms dirichlet(refcat = 2.5)", brms::dirichlet(refcat = 2.5)$refcat)
show("brms dirichlet(refcat = TRUE)", brms::dirichlet(refcat = TRUE)$refcat)
rd <- tools::Rd_db("frmtmb", lib.loc = "C:/Users/adf44/source/r/wt-famsimplex-lib")
txt <- utils::capture.output(tools::Rd2txt(rd[["confint.frmtmb_fit.Rd"]]))
i <- grep("lncor|thetar|chart", txt)
cat("\n--- ?confint lines naming the charts\n")
cat(txt[sort(unique(c(i, i + 1, i - 1)))], sep = "\n")
txt2 <- utils::capture.output(tools::Rd2txt(rd[["dharma_residuals.Rd"]]))
cat("\n--- ?dharma_residuals lines\n")
cat(txt2[grep("simplex|composition|dirichlet", txt2)], sep = "\n")
