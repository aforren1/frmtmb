# frmtmb's post-fit shapes on lane catmulti's cases, against the brms
# 2.23.0 fits of dev/catmulti-brms-fits.R on the same data. Arg 1: "base"
# or "lane" (which library comes first).
# Output: dev/catmulti-out/shapes-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- "C:/Users/adf44/source/r/wt-catmulti-lib"
base <- "C:/Users/adf44/source/r/rellib-r9"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(lib, base, user) else c(base, user))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "from",
    dirname(find.package("frmtmb.sample")), "\n")
source("dev/catmulti-data.R")
out <- "dev/catmulti-out"
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  if (is.character(x)) return(cat(nm, ":", x, "\n"))
  cat(nm, ": dim", paste(dim(x) %||% length(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) {
    cat("  dimnames[[", k, "]]: ",
        if (is.null(dn[[k]])) "NULL" else
          paste(utils::head(dn[[k]], 6), collapse = ", "), "\n", sep = "")
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
tr <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
# z of frmtmb's estimate against the brms posterior: (ML - mean) / sd
agree <- function(fit, bfit) {
  fb <- brms::fixef(bfit)
  ff <- fixef(fit)
  nm <- rownames(fb)
  cat("names identical:", identical(rownames(ff), nm), "\n")
  z <- (ff[nm, "Estimate"] - fb[nm, "Estimate"]) / fb[nm, "Est.Error"]
  print(round(cbind(frmtmb = ff[nm, "Estimate"], brms = fb[nm, "Estimate"],
                    z = z, se_ratio = ff[nm, "Est.Error"] /
                      fb[nm, "Est.Error"]), 4))
}
post <- function(fit) {
  shape("fitted", tr(fitted(fit)))
  shape("fitted(scale = 'linear')", tr(fitted(fit, scale = "linear")))
  shape("predict", tr(predict(fit)))
}

dc <- catmulti_cat_data(11)
sec("cat-refb")
f1 <- tr(frm(bf(f ~ x), family = categorical(refcat = "b"), data = dc))
if (!is.character(f1)) {
  agree(f1, readRDS(file.path(out, "brms-cat-refb.rds")))
  post(f1)
  print(fitted(f1)[1:2, , ])
}
sec("brms::categorical(refcat = 'b')")
f1b <- tr(frm(bf(f ~ x), family = brms::categorical(refcat = "b"),
              data = dc))
if (!is.character(f1b)) print(rownames(fixef(f1b)))
sec("cat-refna")
w <- character(0)
f2 <- withCallingHandlers(
  tr(frm(bf(f ~ x), family = categorical(refcat = NA), data = dc)),
  warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
cat("warnings:", length(w), "\n"); cat(substr(w, 1, 300), sep = "\n")
if (!is.character(f2)) {
  print(fixef(f2))
  f0 <- frm(bf(f ~ x), family = categorical(), data = dc)
  cat("logLik NA - logLik ref a:",
      format(as.numeric(logLik(f2) - logLik(f0)), digits = 4), "\n")
  cat("max |fitted NA - fitted ref a|:",
      format(max(abs(fitted(f2)[, 1, ] - fitted(f0)[, 1, ])), digits = 4),
      "\n")
  post(f2)
}
dm <- catmulti_mult_data(12)
sec("mult-bare")
f3 <- tr(frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dm))
if (!is.character(f3)) {
  agree(f3, readRDS(file.path(out, "brms-mult-bare.rds")))
  post(f3)
  print(tr(fitted(f3)[1:2, , ]))
  shape("predict", tr(predict(f3)))
  nd <- data.frame(x = c(-1, 0, 1), N = c(5L, 10L, 20L))
  print(tr(fitted(f3, newdata = nd)))
  print(tr(fitted(f3, newdata = data.frame(x = c(-1, 0, 1)))))
}
sec("mult-refy2")
f4 <- tr(frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = "y2"),
             data = dm))
if (!is.character(f4)) agree(f4, readRDS(file.path(out, "brms-mult-refy2.rds")))
sec("default_prior")
print(tr(default_prior(bf(f ~ x), data = dc,
                       family = categorical(refcat = "b"))))
print(tr(default_prior(bf(Yc | trials(N) ~ x), data = dm,
                       family = multinomial())))
sec("frm_sample draws: categorical refcat b")
set.seed(5)
ds <- tr(frm_sample(f1, chains = 1, iter = 300, warmup = 150,
                    refresh = 0, seed = 5))
if (!is.character(ds)) {
  shape("posterior_linpred", tr(posterior_linpred(ds)))
  shape("posterior_linpred(dpar = 'mua')",
        tr(posterior_linpred(ds, dpar = "mua")))
  shape("posterior_linpred(transform = TRUE)",
        tr(posterior_linpred(ds, transform = TRUE)))
  shape("posterior_epred", tr(posterior_epred(ds)))
}
sec("frm_sample draws: multinomial")
ds3 <- tr(frm_sample(f3, chains = 1, iter = 300, warmup = 150,
                     refresh = 0, seed = 5))
if (!is.character(ds3)) {
  shape("posterior_linpred", tr(posterior_linpred(ds3)))
  shape("posterior_epred", tr(posterior_epred(ds3)))
  shape("posterior_predict", tr(posterior_predict(ds3)))
}
