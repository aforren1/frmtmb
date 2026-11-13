# Lane ordstart, item 2: which part of the degenerate rule decides the
# ranking. Over the 40-start runs (dev/ordstart-starts.R, nstart 40), per
# seed: the best candidate overall, the best that no part of the warning's
# rule flags ("rule-proper"), and the best with no step-function
# component ("step-free": sharpen <= -0.1 in every component), each with
# its logLik and kind; and the cover histogram (smallest, over a
# component's categories, of the largest probability any row gives).
#   Rscript dev/ordstart-classes.R file.tsv [file.tsv ...]
args <- commandArgs(trailingOnly = TRUE)
X <- do.call(rbind, lapply(args, utils::read.delim))
X <- X[!is.na(X$ll), ]
X$step <- pmax(X$sharpen1, X$sharpen2, na.rm = TRUE) > -0.1
X$runoff <- pmax(X$near1, X$near2, na.rm = TRUE) > 50
X$coll <- X$collapsed1 | X$collapsed2
X$cover <- pmin(X$cover1, X$cover2)
X$kind <- ifelse(X$step, "step", ifelse(X$runoff | X$coll,
                 paste0(ifelse(X$runoff, "runoff", ""),
                        ifelse(X$coll, "collapsed", "")), "proper"))
cat("kinds over", nrow(X), "runs:\n")
print(table(X$kind, X$family))
cat("\ncover (log10) by kind:\n")
h <- cut(log10(pmax(X$cover, 1e-320)),
         c(-400, -100, -30, -16, -12, -10, -8, -6, -2, 0))
print(table(h, X$kind))
out <- NULL
for (key in unique(paste(X$family, X$seed))) {
  S <- X[paste(X$family, X$seed) == key, ]
  d0 <- S[S$run == 0, ]
  b <- S[which.max(S$ll), ]
  R <- S[S$kind == "proper", ]
  P <- S[!S$step, ]
  out <- rbind(out, data.frame(
    family = S$family[1], seed = S$seed[1], ll0 = d0$ll, kind0 = d0$kind,
    best = b$ll, kind_best = b$kind,
    rule_proper = if (nrow(R)) max(R$ll) else NA,
    step_free = if (nrow(P)) max(P$ll) else NA,
    kind_step_free = if (nrow(P)) P$kind[which.max(P$ll)] else NA,
    n_step = sum(S$step), n_runoff_only = sum(S$kind == "runoff"),
    n_coll_only = sum(S$kind %in% c("collapsed", "runoffcollapsed"))))
}
out$sf_minus_rp <- out$step_free - out$rule_proper
op <- options(width = 200)
print(out, digits = 7, row.names = FALSE)
d <- out$sf_minus_rp
cat("\nseeds", nrow(out), "; step-free best above rule-proper best by > 1e-3:",
    sum(d > 1e-3, na.rm = TRUE), "; equal within 1e-3:",
    sum(abs(d) <= 1e-3, na.rm = TRUE), "; no rule-proper candidate:",
    sum(is.na(out$rule_proper)), "; no step-free candidate:",
    sum(is.na(out$step_free)), "\n")
