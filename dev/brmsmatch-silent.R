## The positional calls this package ANSWERS with something other than
## what brms answers, with nothing said. The formals audit in
## dev/brmsmatch-beyond.R cannot see these: it truncates each side at
## its own `...`, so a method whose arguments run out first scores as
## agreeing and the positional argument lands in `...` and is dropped.
##
## Reproduced on this lane's own cached draws rather than copied from
## dev/reviews/20260915-brmsmatch.md, which found them on its own.
##
##   Rscript dev/brmsmatch-silent.R
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))

ds <- readRDS("dev/stan-cache/brmsmatch-draws.rds")
cat("draws ", nrow(ds$draws), " x ", ncol(ds$draws), "\n\n", sep = "")

shape <- function(v) {
  if (inherits(v, "try-error")) {
    return(paste0("ERROR: ", sub("\n.*$", "",
                                 sub("^Error[^:]*: ", "",
                                     as.character(v)))))
  }
  d <- dim(v)
  paste0(class(v)[[1L]], " ",
         if (is.null(d)) paste0("len=", length(v)) else
           paste(d, collapse = "x"))
}
ev <- function(e) try(suppressWarnings(suppressMessages(e)),
                      silent = TRUE)

## Each row: the positional call, what brms's method would do with that
## slot, and whether the answer here is INDISTINGUISHABLE from the same
## call with no second argument at all.
## narrow enough to paste into dev/brmsmatch-findings.md, where house
## style is 80 columns
cat(sprintf("%-31s %-21s %s\n", "positional call", "what you get",
            "vs no 2nd arg"))
rows <- list(
  list("as.matrix(ds, \"x\")", quote(as.matrix(ds, "x")),
       quote(as.matrix(ds)), "brms: 1 column"),
  list("as.array(ds, \"x\")", quote(as.array(ds, "x")),
       quote(as.array(ds)), "brms: 1 variable"),
  list("as_draws_array(ds, \"x\")", quote(as_draws_array(ds, "x")),
       quote(as_draws_array(ds)), "brms: 1 variable"),
  list("as_draws_df(ds, \"x\")", quote(as_draws_df(ds, "x")),
       quote(as_draws_df(ds)), "brms: 1 variable"),
  list("as_draws_matrix(ds, \"x\")", quote(as_draws_matrix(ds, "x")),
       quote(as_draws_matrix(ds)), "brms: 1 variable"),
  list("fixef(ds, FALSE)", quote(fixef(ds, FALSE)), quote(fixef(ds)),
       "brms: summary = FALSE, RAW DRAWS"),
  list("ranef(ds, FALSE)", quote(ranef(ds, FALSE)), quote(ranef(ds)),
       "brms: summary = FALSE, RAW DRAWS"),
  list("coef(ds, FALSE)", quote(coef(ds, FALSE)), quote(coef(ds)),
       "brms: summary = FALSE, RAW DRAWS"),
  list("VarCorr(ds, NULL, FALSE)", quote(VarCorr(ds, NULL, FALSE)),
       quote(VarCorr(ds)), "brms: summary = FALSE, RAW DRAWS"),
  list("summary(ds, NULL, 0.5)", quote(summary(ds, NULL, 0.5)),
       quote(summary(ds)), "brms: prob = 0.5, 50% interval"),
  list("bayes_R2(ds, NULL, TRUE, TRUE)",
       quote(bayes_R2(ds, NULL, TRUE, TRUE)), quote(bayes_R2(ds)),
       "brms: robust = TRUE, median/MAD")
)
n_silent <- 0L
n_same <- 0L
for (r in rows) {
  a <- ev(eval(r[[2L]]))
  b <- ev(eval(r[[3L]]))
  same <- !inherits(a, "try-error") && !inherits(b, "try-error") &&
    identical(a, b)
  if (!inherits(a, "try-error")) n_silent <- n_silent + 1L
  if (same) n_same <- n_same + 1L
  cat(sprintf("%-31s %-21s %s\n", r[[1L]], shape(a),
              if (same) "IDENTICAL" else "differs"))
}
cat("\npositional calls ANSWERED rather than refused: ", n_silent,
    " of ", length(rows), "\n")
cat("answers indistinguishable from the no-argument call: ", n_same,
    "\n")

## bayes_R2 is one of the six the formals audit DOES see, and it is in
## the severe subclass: it answers rather than erroring.
cat("\n== bayes_R2, the one of the six that answers ==\n")
cat("brms:  bayes_R2.brmsfit(object, resp, summary, robust, probs, ...)\n")
cat("here:  ")
tb <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
          inherits = FALSE)
f <- tryCatch(get("bayes_R2.frmtmb_draws", envir = tb, inherits = FALSE),
              error = function(e) NULL)
if (is.null(f)) {
  for (p in c("rstantools", "brms", "frmtmb")) {
    t <- tryCatch(get(".__S3MethodsTable__.", envir = asNamespace(p),
                      inherits = FALSE), error = function(e) NULL)
    if (is.null(t)) next
    f <- tryCatch(get("bayes_R2.frmtmb_draws", envir = t,
                      inherits = FALSE), error = function(e) NULL)
    if (is.function(f)) break
  }
}
cat("bayes_R2.frmtmb_draws(",
    paste(names(formals(f)), collapse = ", "), ")\n")
r4 <- ev(bayes_R2(ds, NULL, TRUE, TRUE))
if (!inherits(r4, "try-error")) {
  cat("bayes_R2(ds, NULL, TRUE, TRUE) columns: ",
      paste(colnames(r4), collapse = " "), "\n")
  cat("  brms would give a robust median/MAD summary; the fourth\n")
  cat("  positional slot here is `probs`, so TRUE became a quantile\n")
}
cat("DONE\n")
