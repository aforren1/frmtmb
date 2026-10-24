# cox() and royston_parmar() fits on the base build (cox() in core,
# royston_parmar() in frmtmb.spline, rellib-r11) and on the lane's build
# (both in frmtmb.survival), saved so that the two can be compared:
# the objective, the estimates and the standard errors.
#
#   Rscript dev/splitsurv-fit-dump.R <before|after> <out.rds>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
out <- a[2]
LIB <- "C:/Users/adf44/source/r/wt-splitsurv-lib"
REL <- "C:/Users/adf44/source/r/rellib-r11"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "after")) LIB, REL, USR))
if (identical(arm, "after")) {
  suppressMessages(library(frmtmb.survival))
} else {
  suppressMessages(library(frmtmb.spline))
}
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "; cox from",
    environmentName(environment(cox)), "\n")
summ <- function(fit) {
  list(logLik = as.numeric(logLik(fit)), par = fit$opt$par,
       objective = fit$opt$objective, conv = fit$opt$convergence,
       se = tryCatch(sqrt(diag(vcov(fit))), error = function(e) NA),
       fixef = fixef_by_dpar(fit),
       baseline = tryCatch(cox_baseline(fit), error = function(e) NULL),
       mo_search = fit$opt$mo_search)
}
res <- list()
keep <- function(nm, expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(summ(expr), error = function(e) {
    structure(conditionMessage(e), class = "dump_error")
  }), warning = function(z) {
    w <<- c(w, conditionMessage(z))
    invokeRestart("muffleWarning")
  })
  if (is.list(v)) v$warnings <- w
  res[[nm]] <<- v
}
cox_data <- function(seed, n = 300) {
  set.seed(seed)
  x <- rnorm(n)
  g <- factor(sample(1:20, n, replace = TRUE))
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x + rnorm(20, 0, 0.5)[g]))
  ct <- stats::rexp(n, 0.2)
  data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
             g = g, z = rnorm(n))
}
d1 <- cox_data(162)
keep("cox_cens", frm(bf(time | cens(cens) ~ x), family = cox(), data = d1))
keep("cox_face40", frm(bf(time | cens(cens) ~ x), family = cox(),
                       data = cox_data(40)))
keep("cox_frailty", frm(bf(time | cens(cens) ~ x + (1 | g)),
                        family = cox(), data = d1))
keep("cox_df6_identity", frm(bf(time | cens(cens) ~ x),
                             family = cox(link = "identity", df = 6),
                             data = transform(d1, x = abs(x))))
keep("cox_trunc", frm(bf(time | trunc(lb = 1e-4) ~ x), family = cox(),
                      data = d1))
keep("cox_prior", frm(bf(time | cens(cens) ~ x), family = cox(), data = d1,
                      prior = set_prior("dirichlet(2)", class = "sbhaz")))
keep("cox_string", frm(bf(time | cens(cens) ~ x), family = "cox", data = d1))
keep("cox_brmsfamily", frm(bf(time | cens(cens) ~ x),
                           family = brmsfamily("cox"), data = d1))
gp <- get_prior(bf(time | cens(cens) ~ x), data = d1, family = cox())
res$cox_get_prior <- as.data.frame(gp)
e <- new.env()
utils::data("bc", package = "flexsurv", envir = e)
bc <- e$bc
bc$censored <- 1 - bc$censrec
keep("rp_bc", frm(bf(recyrs | cens(censored) ~ group),
                  family = royston_parmar(df = 2), data = bc))
keep("rp_bc_odds", frm(bf(recyrs | cens(censored) ~ group),
                       family = royston_parmar(df = 2, scale = "odds"),
                       data = bc))
saveRDS(res, out)
for (nm in names(res)) {
  v <- res[[nm]]
  cat(sprintf("%-18s %s\n", nm, if (inherits(v, "dump_error")) {
    paste("ERROR:", substr(unclass(v), 1, 80))
  } else if (is.data.frame(v)) "table" else {
    sprintf("logLik %.10f, %d warnings", v$logLik, length(v$warnings))
  }))
}
