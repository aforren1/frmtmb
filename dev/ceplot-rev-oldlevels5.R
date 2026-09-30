# Reviewer re-check (lane ceplot, punch round 1): inside ONE predict()
# call on a multivariate model whose responses both read (1 | g), which
# seen level does each response's unseen level read? The choice is
# recorded by wrapping predict_new_level_spec() in the namespace.
# Data seed 3; calls at seeds 1 to 12.
#   Rscript dev/ceplot-rev-oldlevels5.R > dev/ceplot-rev-log/p1/oldlevels5.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(3)
ng <- 12L
d <- data.frame(x = rnorm(360), g = factor(rep(seq_len(ng), 30)))
u <- rnorm(ng, 0, 2)
d$y <- rnorm(360, 1 + 0.5 * d$x + u[d$g])
d$y2 <- rnorm(360, -1 + 0.3 * d$x + u[d$g] + rnorm(ng)[d$g])
fv <- frm(bf(y ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)), family = gaussian(),
          data = d)
ns <- asNamespace("frmtmb")
orig <- get("predict_new_level_spec", ns)
rec <- new.env()
wrapped <- function(...) {
  out <- orig(...)
  rec$picks <- c(rec$picks, paste(unlist(attr(out, "old_pick")),
                                  collapse = ","))
  out
}
unlockBinding("predict_new_level_spec", ns)
assign("predict_new_level_spec", wrapped, envir = ns)
nd <- data.frame(x = 0, g = factor("n1"))
res <- vapply(1:12, function(s) {
  rec$picks <- character(0)
  set.seed(s)
  invisible(predict(fv, newdata = nd, allow_new_levels = TRUE,
                    sample_new_levels = "old_levels", ndraws = 5))
  paste(rec$picks, collapse = " / ")
}, "")
cat("one predict() call per seed, the choice recorded per response:\n")
for (s in 1:12) cat("  seed", s, ":", res[s], "\n")
