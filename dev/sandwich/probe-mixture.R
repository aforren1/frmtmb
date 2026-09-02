# Probe: does the group-level (latent-class) mixture branch factor over
# clusters? Its per-row weights sit INSIDE the group sum, so a
# masked-out group contributes log(sum_k pi_k) = 0 rather than dropping
# out of the tape - which is the same exact zero the random-effect
# prior gives, but by a different route, so it is worth checking.
#
# Run: Rscript dev/sandwich/probe-mixture.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
con <- file("dev/sandwich/probe-mixture.txt", "w")
say <- function(...) writeLines(paste0(...), con)

set.seed(21)
G <- 30
m <- 5
gv <- factor(rep(seq_len(G), each = m))
cls <- rbinom(G, 1, 0.4)
y <- rnorm(G * m, c(0, 3)[cls[gv] + 1L], 1)
dd <- data.frame(y = y, g = gv)
fit <- try(frm(bf(y ~ 1) +
                 mixture(gaussian(), gaussian(), groups = ~ g),
               data = dd), silent = TRUE)
if (inherits(fit, "try-error")) {
  say("mixture fit failed: ", conditionMessage(attr(fit, "condition")))
} else {
  say("mixture fit ok, mix_g present: ", !is.null(fit$frame$mix_g))
  S <- try(cluster_scores(fit, ~ g), silent = TRUE)
  if (inherits(S, "try-error")) {
    say("cluster_scores failed: ", conditionMessage(attr(S, "condition")))
  } else {
    g0 <- drop(fit$obj$gr(fit$opt$par))
    say("G rows: ", nrow(S))
    say("max |colSums(S) + obj$gr| = ",
        format(max(abs(colSums(S) + g0)), digits = 3))
    say("max |obj$gr| = ", format(max(abs(g0)), digits = 3))
  }
  # a coarser cluster (pairs of groups) must also work
  cl2 <- factor(((as.integer(gv) - 1L) %/% 2L) + 1L)
  S2 <- try(cluster_scores(fit, cl2), silent = TRUE)
  say("coarser cluster: ", if (inherits(S2, "try-error"))
    conditionMessage(attr(S2, "condition")) else
      paste("ok,", nrow(S2), "rows, max err",
            format(max(abs(colSums(S2) +
                             drop(fit$obj$gr(fit$opt$par)))),
                   digits = 3)))
  # a cluster that SPLITS a mixture group must be refused
  cl3 <- factor(seq_len(G * m))
  r3 <- try(cluster_scores(fit, cl3), silent = TRUE)
  say("splitting cluster: ", if (inherits(r3, "try-error"))
    conditionMessage(attr(r3, "condition")) else "NOT REFUSED (bug)")
}
close(con)
