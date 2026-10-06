# Lane fixes, item 2: why the refit of brmsfit-methods:955 stops at its
# default start, and whether any start fits it.
#   Rscript dev/fixes-u955.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
d <- brms_fixture_data(2)
d$Trt <- as.numeric(as.character(d$Trt))
cat("n =", nrow(d), " mean(count) =", mean(d$count),
    " min(count) =", min(d$count), "\n")
# the formula update(fit2, formula. = bf(count ~ a + b, nl = TRUE))
# builds: the new body with fit2's parameter formulas (brms's
# update.brmsformula(mode = "replace") keeps them)
fo <- bf(count | weights(AgeSD) ~ a + b, a ~ Age + (1 | ID1 | patient),
         b ~ Age + (1 | ID1 | patient), nl = TRUE)
fam <- Gamma("identity")
fit2 <- brms_fixture(2)
up_f <- tryCatch(withCallingHandlers(
  update(fit2, formula. = bf(count ~ a + b, nl = TRUE),
         dry_run = "spec"),
  message = function(m) invokeRestart("muffleMessage")),
  error = function(e) conditionMessage(e))
cat("update(dry_run = 'spec') class:", class(up_f)[1], "\n")

# 1. the objective at the default start
ob <- frm(fo, data = d, family = fam, dry_run = "objective")
tpl <- ob$frame$par_template %||% NULL
st <- ob$obj$par
cat("default start (outer):\n"); print(st)
cat("nll at default start:", ob$obj$fn(st), "\n")
g <- tryCatch(ob$obj$gr(st), error = function(e) conditionMessage(e))
cat("gradient at default start:", format(g), "\n")

# 2. a start with mu > 0: a_Intercept = mean(count) / 2 = b_Intercept
pt <- par_template(fo, data = d, family = fam)
cat("par_template beta names:", names(pt$beta), "\n")
b0 <- pt$beta
b0[] <- 0
b0[grep("Intercept", names(b0))] <- mean(d$count) / 2
fitB <- tryCatch(withCallingHandlers(
  frm(fo, data = d, family = fam, start = list(beta = b0)),
  warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  })
if (!is.null(fitB)) {
  cat("converged:", fitB$opt$convergence, " message:", fitB$opt$message,
      "\n")
  cat("logLik:", format(as.numeric(logLik(fitB)), digits = 12), "\n")
  print(fixef(fitB))
  # the ridge: shift t from b_Intercept into a_Intercept (and the Age
  # slopes) and evaluate the marginal objective; the body reads a + b
  p <- fitB$obj$env$last.par.best[-fitB$obj$env$random]
  nm <- names(p)
  bi <- which(nm == "beta")
  cat("beta names in the outer vector:", names(fitB$estimates$beta), "\n")
  base <- fitB$obj$fn(p)
  for (tt in c(-3, -1, 1, 3)) {
    q <- p
    q[bi[1]] <- q[bi[1]] + tt
    q[bi[3]] <- q[bi[3]] - tt
    cat(sprintf("  shift %+d in (a_Int, b_Int): nll %.12f  diff %.3g\n",
                tt, fitB$obj$fn(q), fitB$obj$fn(q) - base))
  }
  H <- optimHess(p, fitB$obj$fn, fitB$obj$gr)
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  cat("outer Hessian eigenvalues:", format(ev, digits = 3), "\n")
  cat("ratio smallest / largest:", format(min(abs(ev)) / max(abs(ev)),
                                         digits = 3), "\n")
}

# 3. brms's own fit2 priors, normal(2, 2) on a and normal(0, 3) on b,
# which brms's update() carries over: a penalty here, and a start
pr <- c(set_prior("normal(2, 2)", nlpar = "a"),
        set_prior("normal(0, 3)", nlpar = "b"))
fitC <- tryCatch(withCallingHandlers(
  frm(fo, data = d, family = fam, prior = pr),
  warning = function(w) {
    cat("  [warning]", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("  [message]", conditionMessage(m)); invokeRestart("muffleMessage")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  })
if (!is.null(fitC)) {
  cat("with brms priors: converged:", fitC$opt$convergence, "\n")
  print(fixef(fitC))
}
