# Reviewer: brms 2.23.0's answers on the lane's cached fits, for the
# routes the lane refuses or reshapes. No compile needed.
# Output: dev/famsimplex-rev-brms.log
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
shape <- function(lab, expr) {
  cat("\n---", lab, "\n")
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  if (length(w)) cat("WARNINGS:", unique(w), sep = "\n  ")
  if (is.character(r) && length(r) == 1) {
    cat(r, "\n")
  } else if (inherits(r, "brms_conditional_effects")) {
    d <- r[[1]]
    cat("names:", names(d), "\n")
    cat("nrow:", nrow(d), " cats__:", levels(d$cats__), "\n")
    print(utils::head(d[, intersect(c("x", "cats__", "estimate__", "se__",
                                      "lower__", "upper__"), names(d))], 3))
  } else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
  }
  invisible(r)
}
set.seed(1)
for (nm in c("dir", "ln", "dir2", "dirre")) {
  bfit <- readRDS(file.path(wt, paste0("famsimplex-brmsfit-", nm, ".rds")))
  cat("\n==========", nm, "==========\n")
  nd <- bfit$data[c(5, 1, 3), setdiff(names(bfit$data), "Y")]
  nd$Y <- bfit$data$Y[c(5, 1, 3), ]
  shape("fitted newdata", fitted(bfit, newdata = nd, ndraws = 200))
  shape("fitted newdata dpar mub", fitted(bfit, newdata = nd, ndraws = 200,
                                          dpar = "mub"))
  shape("posterior_epred dpar mub",
        posterior_epred(bfit, newdata = nd, ndraws = 20, dpar = "mub"))
  shape("predict newdata", predict(bfit, newdata = nd, ndraws = 200))
  shape("posterior_predict newdata",
        posterior_predict(bfit, newdata = nd, ndraws = 20))
  shape("residuals pearson", residuals(bfit, type = "pearson", ndraws = 20))
  shape("residuals method=posterior_predict",
        residuals(bfit, method = "posterior_predict", ndraws = 20))
  shape("pp_check type stat", pp_check(bfit, type = "stat", ndraws = 5))
  shape("CE x", conditional_effects(bfit, "x"))
  shape("CE x categorical", conditional_effects(bfit, "x",
                                                categorical = TRUE))
  shape("CE x method=posterior_predict",
        conditional_effects(bfit, "x", method = "posterior_predict"))
  shape("CE x method=posterior_predict categorical",
        conditional_effects(bfit, "x", method = "posterior_predict",
                            categorical = TRUE))
  shape("CE x method=posterior_linpred",
        conditional_effects(bfit, "x", method = "posterior_linpred"))
  shape("CE x dpar mub", conditional_effects(bfit, "x", dpar = "mub"))
  if (startsWith(nm, "ln")) {
    shape("CE x dpar sigmab", conditional_effects(bfit, "x",
                                                  dpar = "sigmab"))
  } else if (nm != "dir2") {
    shape("CE x dpar phi", conditional_effects(bfit, "x", dpar = "phi"))
  }
  shape("loo_predict", loo_predict(bfit))
  shape("posterior_linpred transform=TRUE",
        posterior_linpred(bfit, transform = TRUE, ndraws = 5))
  shape("bayes_R2", bayes_R2(bfit))
  shape("prior_summary", {
    ps <- prior_summary(bfit)
    print(ps[, c("prior", "class", "dpar", "lb")])
    "printed"
  })
  cat("variables:", paste(variables(bfit), collapse = " "), "\n")
}
