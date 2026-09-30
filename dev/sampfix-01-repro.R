# Lane sampfix, script 01: reproduce the six filed defects on one build.
#
#   Rscript dev/sampfix-01-repro.R <lane|ref>
#
# lane = C:/Users/adf44/source/r/wt-sampfix-lib first, then rellib-r3;
# ref = rellib-r3 alone. Every block is wrapped so one failure does not
# hide the next.

arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("ARM ", arm, " frmtmb.sample from ",
    dirname(find.package("frmtmb.sample")), "\n", sep = "")

try_show <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) {
    cat("  [", lab, "] ERROR: ", conditionMessage(r), "\n", sep = "")
  } else if (is.numeric(r)) {
    cat("  [", lab, "] numeric ", paste(dim(r) %||% length(r), collapse = "x"),
        ", non-finite ", sum(!is.finite(r)), " of ", length(r), "\n",
        sep = "")
  } else {
    cat("  [", lab, "] OK ", paste(class(r), collapse = "/"), "\n", sep = "")
  }
  if (length(w)) cat("      warnings: ", length(w), " (",
                     length(unique(w)), " distinct): ",
                     unique(w)[1L], "\n", sep = "")
  invisible(r)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(1212L)
dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.7)

cat("\n== 1. laplace draws, predictive methods\n")
for (form in list(bf(y ~ x + ar(t, g) + (1 | g)), bf(y ~ x + (1 | g)))) {
  cat("-- ", deparse(form$formula), "\n", sep = "")
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = gaussian(), data = dd, chains = 1,
               iter = 300, refresh = 0, seed = 3, laplace = TRUE)))
  cat("  columns: ", paste(colnames(ds$draws), collapse = " "), "\n", sep = "")
  set.seed(9)
  try_show("posterior_predict", posterior_predict(ds))
  try_show("posterior_epred", posterior_epred(ds))
  try_show("posterior_linpred", posterior_linpred(ds))
  try_show("posterior_epred re_formula=NA",
           posterior_epred(ds, re_formula = NA))
  try_show("log_lik", log_lik(ds))
}

cat("\n== 2. laplace = TRUE with no random effect\n")
try_show("frm_sample(y ~ x, laplace = TRUE)",
         frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 1,
                    iter = 300, refresh = 0, seed = 3, laplace = TRUE))

cat("\n== 4. stanfit = NULL draws object\n")
fit0 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
dsr <- suppressWarnings(suppressMessages(
  frm_sample(fit0, chains = 2, iter = 200, refresh = 0, seed = 3)))
dn <- dsr
dn$stanfit <- NULL
try_show("nchains", nchains(dn))
try_show("ndraws", ndraws(dn))
try_show("VarCorr", VarCorr(dn))
try_show("as_draws_array", as_draws_array(dn))
try_show("summary", summary(dn))

cat("\n== 5. variables() of ordinal and cs() draws\n")
set.seed(405L)
n <- 300L
do <- data.frame(x = rnorm(n), fc = factor(sample(c("a", "b", "c"), n, TRUE)))
eta <- 0.8 * do$x
do$yo <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.7, Inf), labels = FALSE)
do$yo <- factor(do$yo, ordered = TRUE)
fo <- frm(bf(yo ~ x), family = sratio(), data = do)
dso <- suppressWarnings(suppressMessages(
  frm_sample(fo, chains = 1, iter = 200, refresh = 0, seed = 3)))
cat("  variables(fit): ", paste(variables(fo), collapse = " "), "\n", sep = "")
cat("  variables(ds):  ", paste(variables(dso), collapse = " "), "\n", sep = "")
fc <- frm(bf(yo ~ x + cs(fc)), family = sratio(), data = do)
dsc <- suppressWarnings(suppressMessages(
  frm_sample(fc, chains = 1, iter = 200, refresh = 0, seed = 3)))
cat("  variables(fit): ", paste(variables(fc), collapse = " "), "\n", sep = "")
cat("  variables(ds):  ", paste(variables(dsc), collapse = " "), "\n", sep = "")

cat("\n== 6. one-parameter model\n")
fx <- frm(bf(y ~ 0 + x), family = poisson(),
          data = transform(dd, y = rpois(nrow(dd), 2)))
try_show("frm_sample(y ~ 0 + x) poisson",
         frm_sample(fx, chains = 1, iter = 200, refresh = 0, seed = 3))
cat("\nDONE\n")
