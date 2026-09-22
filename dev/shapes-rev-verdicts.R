# Reviewer, priority 6: do the 40 dropped verdict rows pass for the
# RIGHT reason? A shape assertion can be satisfied by a matrix of NA,
# and a row-equality assertion is vacuous when both rows are NA. So
# print the OBJECT each dropped assertion reads, not just the verdict.
#
#   Rscript dev/shapes-rev-verdicts.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
setwd(TREE)
# the fixtures are the ported tier's own, read from its helper
source("tests/testthat/helper-brms-suite.R")

fitk <- function(k) suppressWarnings(brms_fixture(k))
show <- function(id, lab, e) {
  v <- tryCatch(suppressWarnings(e),
                error = function(c) paste0("ERROR: ",
                                           conditionMessage(c)))
  cat("\n-- ", id, "  ", lab, "\n", sep = "")
  if (is.character(v) && length(v) == 1 && startsWith(v, "ERROR")) {
    cat("   ", v, "\n"); return(invisible(NULL))
  }
  if (is.array(v) || is.matrix(v)) {
    cat("   dim ", paste(dim(v), collapse = " x "), "  finite cells ",
        sum(is.finite(as.numeric(v))), " of ", length(v), "\n", sep = "")
    if (length(dim(v)) == 2L) print(utils::head(round(v, 5), 3))
    else print(round(v[seq_len(min(2, dim(v)[1])), , , drop = FALSE], 5))
  } else {
    print(v)
  }
  invisible(v)
}

cat("############ fixture 1 ############\n")
f1 <- fitk(1)
cat("convergence: ", f1$opt$convergence, "  message: ",
    paste(f1$opt$message, collapse = " "), "\n")
fi <- show(":292 :293", "fitted(fit1)", fitted(f1))
show(":364", "rownames(fixef(fit1))", rownames(fixef(f1)))
show(":364b", "fixef(fit1) finite Est.Error",
     sum(is.finite(fixef(f1)[, "Est.Error"])))
show(":584", "ngrps(fit1)", ngrps(f1))
show(":784", "grep 'Multilevel Hyperparameters:' in print(fit1)",
     any(grepl("Multilevel Hyperparameters:",
               utils::capture.output(print(f1)))))
s1 <- suppressWarnings(summary(f1))
show(":887", "is.data.frame(summary1$fixed)", is.data.frame(s1$fixed))
show(":888", "rownames(summary1$fixed)", rownames(s1$fixed))
show(":894", "rownames(summary1$random$visit)",
     rownames(s1$random$visit))
show(":896", "grep 'Regression Coefficients:' in print(summary1)",
     any(grepl("Regression Coefficients:",
               utils::capture.output(print(s1)))))
show(":897", "grep 'Priors:' in print(summary1, priors = TRUE)",
     any(grepl("Priors:", utils::capture.output(
       print(suppressWarnings(summary(f1, priors = TRUE)))))))
show(":999", "dim(vcov(fit1)) and finite cells", vcov(f1))
show(":1000", "dim(vcov(fit1, cor = TRUE))",
     vcov(f1, correlation = TRUE))
show(":726 :727", "predict(fit1)", predict(f1, ndraws = 200))
show(":729", "predict(fit1, probs = c(0.1,0.5,0.9))",
     predict(f1, ndraws = 200, probs = c(0.1, 0.5, 0.9)))
show(":825", "residuals(fit1, probs = 0.5)  [expects n x 3]",
     residuals(f1, probs = 0.5))
show(":326", "fitted(fit1, dpar = 'sigma')", fitted(f1, dpar = "sigma"))

cat("\n############ fixture 2 ############\n")
f2 <- fitk(2)
cat("convergence: ", f2$opt$convergence, "\n")
show(":334 :339 :342", "fitted(fit2)", fitted(f2))
nd <- data.frame(Age = c(0, 0), AgeSD = c(1, 1), Trt = c(1, 1),
                 patient = c(1, 1), count = c(20, 20))
show(":337 :340", "fitted(fit2, newdata = 2 identical rows)",
     fitted(f2, newdata = nd))
show(":340b", "rows 1 and 2 equal AND finite?", {
  v <- fitted(f2, newdata = nd)
  c(equal = isTRUE(all.equal(v[1, ], v[2, ])),
    finite = all(is.finite(v[1, ])))
})
show(":585", "ngrps(fit2)", ngrps(f2))
show(":835", "residuals(fit2)", residuals(f2))
show(":750 :753", "predict(fit2)", predict(f2, ndraws = 200))

cat("\n############ fixture 4 (ordinal) ############\n")
f4 <- fitk(4)
show(":348", "fitted(fit4)", fitted(f4))
show(":761 :762", "predict(fit4)", predict(f4, ndraws = 200))
show(":369", "rownames(fixef(fit4))", rownames(fixef(f4)))

cat("\n############ fixture 5 (mixture) ############\n")
f5 <- fitk(5)
show(":355", "fitted(fit5)", fitted(f5))
show(":767", "predict(fit5)", predict(f5, ndraws = 200))
