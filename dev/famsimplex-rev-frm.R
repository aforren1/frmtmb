# Reviewer: the lane build's answers on the same data and calls as
# dev/famsimplex-rev-brms.R. Output: dev/famsimplex-rev-frm.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
shape <- function(lab, expr) {
  cat("\n---", lab, "\n")
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  if (length(w)) cat("WARNINGS:", unique(w), sep = "\n  ")
  if (is.character(r) && length(r) == 1) {
    cat(r, "\n")
  } else if (inherits(r, "frmtmb_conditional_effects") ||
               (is.list(r) && is.data.frame(r[[1]]))) {
    d <- r[[1]]
    cat("names:", names(d), "\n")
    cat("nrow:", nrow(d), " cats__:", levels(d$cats__), "\n")
    print(utils::head(d[, intersect(c("x", "cats__", "estimate__", "se__",
                                      "lower__", "upper__"), names(d))], 3))
  } else {
    cat("class:", class(r), " dim:", dim(r), "\n")
    if (!is.null(dimnames(r))) str(dimnames(r))
    if (is.null(dim(r))) print(utils::head(r))
  }
  invisible(r)
}
forms <- list(dir = bf(Y ~ x, phi ~ z), dirre = bf(Y ~ x + (1 | ID | g)),
              ln = bf(Y ~ x, sigmab ~ z), dir2 = bf(Y ~ x))
fams <- list(dir = dirichlet(), dirre = dirichlet(), ln = logistic_normal(),
             dir2 = brmsfamily("dirichlet2"))
set.seed(1)
for (nm in c("dir", "ln", "dir2", "dirre")) {
  dat <- readRDS(file.path(wt, paste0("famsimplex-brms-", nm, ".rds")))$data
  fit <- frm(forms[[nm]], family = fams[[nm]], data = dat)
  cat("\n==========", nm, "==========\n")
  nd <- dat[c(5, 1, 3), setdiff(names(dat), "Y")]
  nd$Y <- dat$Y[c(5, 1, 3), ]
  shape("fitted newdata", fitted(fit, newdata = nd))
  shape("fitted newdata dpar mub", fitted(fit, newdata = nd, dpar = "mub"))
  shape("predict newdata", predict(fit, newdata = nd, ndraws = 200))
  shape("simulate rowsums", {
    s <- simulate(fit, nsim = 3, seed = 2)
    range(vapply(s, function(m) max(abs(rowSums(m) - 1)), 0))
  })
  shape("residuals pearson", residuals(fit, type = "pearson"))
  shape("pp_check type stat", pp_check(fit, type = "stat"))
  shape("CE x", conditional_effects(fit, "x"))
  shape("CE x categorical", conditional_effects(fit, "x",
                                                categorical = TRUE))
  shape("CE x categorical FALSE", conditional_effects(fit, "x",
                                                      categorical = FALSE))
  shape("CE x method=posterior_predict",
        conditional_effects(fit, "x", method = "posterior_predict"))
  shape("CE x method=posterior_linpred",
        conditional_effects(fit, "x", method = "posterior_linpred"))
  shape("CE x dpar mub", conditional_effects(fit, "x", dpar = "mub"))
  if (startsWith(nm, "ln")) {
    shape("CE x dpar sigmab", conditional_effects(fit, "x", dpar = "sigmab"))
  } else if (nm != "dir2") {
    shape("CE x dpar phi", conditional_effects(fit, "x", dpar = "phi"))
  }
  shape("frm_linpred type response", frm_linpred(fit, type = "response"))
  shape("bayes_R2", bayes_R2(fit))
  shape("confint", confint(fit))
  shape("summary", {
    print(summary(fit))
    "printed"
  })
  shape("variables", variables(fit))
  shape("emmeans dpar mub", {
    print(emmeans::emmeans(fit, ~ 1, dpar = "mub"))
    "printed"
  })
}
