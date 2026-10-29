# Reviewer: does the instrument have power? Each mutation changes the
# model frmtmb fits (frm side only) while brms keeps the lane's spelling.
# A sound instrument scores every mutant "defect"; an "equivalent" is a
# blind spot. One mutant per R process: Rscript mutate.R <name>.
args <- commandArgs(trailingOnly = TRUE)
which <- args[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb); library(testthat)
  requireNamespace("brms"); requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
# frm side reads spec$prior_frm when the mutant sets it
orig_fit <- E$eq_frm_fit
E$eq_frm_fit <- function(spec, ctrl = NULL) {
  if (!is.null(spec$prior_frm)) spec$prior <- spec$prior_frm
  orig_fit(spec, ctrl)
}
environment(E$eq_frm_fit) <- E
smp <- read.delim("dev/equiv-data/sample.tsv")
M <- list(
  prior_scale = list(key = "gaussian x prior", f = function(s) {
    s$prior_frm <- sub("normal(0, 0.5)", "normal(0, 0.55)", s$prior,
                       fixed = TRUE); s }),
  nlprior = list(key = "prior x nl", f = function(s) {
    s$prior_frm <- sub("normal(0, 1)", "normal(0, 1.5)", s$prior,
                       fixed = TRUE); s }),
  arma_cov = list(key = "me() x arma()", f = function(s) {
    s$bf_frm <- sub("gr = g)", "gr = g, cov = TRUE)", s$bf, fixed = TRUE)
    s }),
  me_sd = list(key = "us x me()", f = function(s) {
    s$data$sx2 <- 0.33
    s$bf_frm <- sub("me(xo, sx)", "me(xo, sx2)", s$bf, fixed = TRUE); s }),
  nl_body = list(key = "zero_inflated_beta x nl", f = function(s) {
    s$bf_frm <- sub("b2 * x", "b2 * x * 1.05", s$bf, fixed = TRUE); s }),
  mo_link = list(key = "poisson x mo()", f = function(s) {
    s$family_frm <- "poisson(link = 'sqrt')"; s }),
  mm_w = list(key = "rescor x mm()", f = function(s) {
    s$data$w1 <- 0.6; s$data$w2 <- 0.4
    s$bf_frm <- sub("mm(g1, g2)", "mm(g1, g2, weights = cbind(w1, w2))",
                    s$bf, fixed = TRUE); s }),
  cens_rows = list(key = "gaussian x cens()", f = function(s) {
    s$data$c2 <- s$data$c
    s$data$c2[which(s$data$c2 == 1L)[1:3]] <- 0L
    s$bf_frm <- sub("cens(c)", "cens(c2)", s$bf, fixed = TRUE); s }),
  thres_prior = list(key = "thres() x prior", f = function(s) {
    s$prior_frm <- sub("normal(0, 2)", "normal(0, 2.5)", s$prior,
                       fixed = TRUE); s }),
  probit_curve = list(key = "cumulative x frm_curve", f = function(s) {
    s$family_frm <- "cumulative(link = 'probit')"; s }),
  sim_sigma = list(key = "gaussian x simulate", f = function(s) {
    # same model; the simulate method is given one wrong row in 200
    s }),
  cs_drop = list(key = "cs_pred() x prior", f = function(s) {
    s$bf_frm <- "bf(y ~ x + z)"; s }),
  # L4 power on a penalized nl model (p042): the smooth replaced by a
  # line, and the smooth's basis cut to k = 4
  nls_line = list(key = "s() x nl", f = function(s) {
    s$bf_frm <- sub("b1 ~ s(z)", "b1 ~ z", s$bf, fixed = TRUE); s }),
  nls_k4 = list(key = "s() x nl", f = function(s) {
    s$bf_frm <- sub("b1 ~ s(z)", "b1 ~ s(z, k = 4)", s$bf, fixed = TRUE)
    s }),
  nls_quad = list(key = "smooth x nl", f = function(s) {
    s$bf_frm <- sub("b1 ~ s(z)", "b1 ~ z + I(z^2)", s$bf, fixed = TRUE)
    s }),
  # D2: is frmtmb's named remedy (class "sd", group "s(x)") brms's
  # class "sds" prior? brms keeps "sds"; frm() reads the remedy.
  d2_remedy = list(key = "s() x prior", f = function(s) {
    s$prior_frm <- paste0("set_prior('student_t(3, 0, 1)', class = 'sd', ",
                          "group = 's(x)')")
    s }))
if (which == "d2_remedy") {
  # sds is a positive block (log chart), so its Jacobian is added as for sd
  opb <- E$eq_prior_blocks
  E$eq_prior_blocks <- function(spec, stan_names) {
    out <- opb(spec, stan_names)
    p <- E$eq_eval(spec$prior, "brms")
    if ("sds" %in% p$class) out <- c(out, grep("^sds_", stan_names,
                                               value = TRUE))
    unique(out)
  }
}
m <- M[[which]]
base <- E$eq_specs[[m$key]]
E$eq_specs[[m$key]] <- function() m$f(base())
row <- smp[paste(smp$feature_a, smp$feature_b, sep = " x ") == m$key, ]
rec <- E$eq_measure(m$key, row)
flat <- E$eq_flatten(rec, row)
cat(sprintf("MUTANT %s (%s %s): verdict=%s level=%s l1_const=%s l1_grad=%s l1_extra=%s l2=%s %s | frm_err=%s | l1_err=%s\n",
            which, row$id, m$key, flat$verdict, flat$level,
            signif(flat$l1_const, 3), signif(flat$l1_grad, 3), flat$l1_extra,
            flat$l2_ok, substr(flat$l2_what, 1, 80),
            flat$l4_z, signif(flat$l4_rhat, 3), substr(flat$frm_error, 1, 100), substr(flat$l1_error, 1, 100)))
