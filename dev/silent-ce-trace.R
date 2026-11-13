# Item 720: which call refuses the vreal() payload in each
# conditional_effects() path, on the base build.
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) -
      0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vreal1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
d <- data.frame(x = rnorm(n))
d$y <- (1 + 0.5 * d$x) * rep_len(c(1, 3), n) + rnorm(n, 0, 0.3)
f_expr <- frm(bf(y | vreal(rep_len(c(1, 3), 120)) ~ x), family = scaled,
              data = d)
trace("aterms_for_newdata", where = asNamespace("frmtmb"), print = FALSE,
      tracer = quote({
        cs <- vapply(sys.calls(), function(cl) deparse1(cl[[1L]])[1L], "")
        cat("  aterms_for_newdata from:",
            paste(utils::tail(cs, 6), collapse = " < "), "\n")
      }))
trace("ce_aterms", where = asNamespace("frmtmb"), print = FALSE,
      tracer = quote(cat("  ce_aterms reached\n")))
# the seam frmtmb.sample's conditional_effects(method = "predict") draws
# through, one call per posterior draw: it reaches ce_aterms() and no
# aterms_for_newdata()
nd <- data.frame(x = seq(-2, 2, length.out = 100), y = mean(d$y))
r <- withCallingHandlers(
  tryCatch(frmtmb:::ce_sim_one(f_expr, nd, "y"),
           error = function(e) paste("ERROR:", conditionMessage(e))),
  warning = function(w) {
    cat("  warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
cat("== ce_sim_one on a 100-row grid: class", class(r), "length", length(r),
    " NA", sum(is.na(r)), " head", paste(format(head(r, 4)), collapse = " "),
    "\n")
cat("   ce_aterms names:",
    paste(names(frmtmb:::ce_aterms(f_expr$spec$responses$y, nd, 100L)),
          collapse = ","), "(empty means dropped)\n")
for (m in c("posterior_epred", "predict")) {
  cat("== method", m, "\n")
  withCallingHandlers(
    tryCatch(conditional_effects(f_expr, effects = "x", method = m),
             error = function(e) cat("ERROR:", conditionMessage(e), "\n")),
    error = function(e) {
      cs <- vapply(sys.calls(), function(cl) deparse1(cl[[1L]]), "")
      cat("  stack:", paste(utils::tail(cs[!grepl("^(withCallingHandlers|tryCatch|doTryCatch|tryCatchList|tryCatchOne|h|signalCondition|stop|frm_stop|cnd_signal)$", cs)], 8), collapse = " < "), "\n")
    })
}
