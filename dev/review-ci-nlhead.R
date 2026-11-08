# Reviewer probe for item 2, the quantity that decides the warning. The
# ubuntu warning fires when the projected gradient exceeds grad_tol (1e-3)
# AND the Newton headroom grad_verdict() computes exceeds grad_tol. The
# gradient alone is a trip-wire. So for each fit, force the headroom
# (grad_tol 1e-12 for the verdict only) and report it against 1e-3.
#
#   Rscript dev/review-ci-nlhead.R <lib1;lib2;...> <base test> <after test>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[[1L]], ";", fixed = TRUE)[[1L]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), "\n")
load_env <- function(path) {
  src <- parse(path)
  keep <- vapply(src, function(e) {
    is.call(e) && identical(e[[1L]], as.name("<-"))
  }, NA)
  env <- new.env(parent = asNamespace("frmtmb"))
  for (e in src[keep]) eval(e, env)
  env
}
eb <- load_env(args[[2L]])
ea <- load_env(args[[3L]])
changed <- c("student", "negbinomial", "negbinomial2", "nbinom1",
             "zero_inflated_negbinomial", "hurdle_negbinomial")
fit1 <- function(env, form, f, d) {
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(form, family = env$nl_fam_of(f), data = d)),
    warning = function(cnd) {
      w <<- c(w, substr(conditionMessage(cnd), 1, 30))
      invokeRestart("muffleWarning")
    })
  v <- frmtmb:::grad_verdict(fit, control = list(grad_tol = 1e-12))
  data.frame(gmax = v$gmax, headroom = v$headroom,
             warn = paste(unique(w), collapse = " | "))
}
rows <- list()
plan <- rbind(data.frame(arm = "before", family = changed),
              data.frame(arm = "after", family = ea$nl_fams_single))
for (i in seq_len(nrow(plan))) {
  f <- plan$family[i]
  env <- if (plan$arm[i] == "before") eb else ea
  z <- env$nl_fam_data()
  d <- data.frame(x = z$x, ntr = 10L)
  d$y <- env$nl_fam_response(f, z)
  L <- env$nl_fam_lhs(f)
  fs <- list(
    lin = bf(stats::as.formula(paste(L, "~ x"))),
    int = bf(stats::as.formula(paste(L, "~ 1"))),
    nl = bf(stats::as.formula(paste(L, "~ b0 + b1 * x")), b0 + b1 ~ 1,
            nl = TRUE))
  for (m in names(fs)) {
    r <- tryCatch(fit1(env, fs[[m]], f, d), error = function(e)
      data.frame(gmax = NA, headroom = NA, warn = conditionMessage(e)))
    rows[[length(rows) + 1L]] <- cbind(arm = plan$arm[i], family = f,
                                       model = m, r)
  }
}
out <- do.call(rbind, rows)
options(width = 200)
out$head_over_tol <- out$headroom / 1e-3
print(format(out, digits = 3), row.names = FALSE)
cat("\nlargest headroom/1e-3 by arm:\n")
print(tapply(out$head_over_tol, out$arm, max, na.rm = TRUE))
a <- out[out$arm == "after", ]
a <- a[order(-a$head_over_tol), ]
cat("\nafter, top 8 by headroom/1e-3:\n")
print(format(utils::head(a, 8), digits = 3), row.names = FALSE)
