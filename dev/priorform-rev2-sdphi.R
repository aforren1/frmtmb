# Reviewer recheck round 1: which density does the phi block's sd get
# when class sd is set class-wide and for dpar phi, in both orders? Read
# off resolve_prior_input() and validate_prior(), on base and lane.
#   Rscript dev/priorform-rev2-sdphi.R ref|lane      seed 20260916
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (mode == "lane") "C:/Users/adf44/source/r/priorform-lib" else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260916)
n <- 200
d <- data.frame(g = factor(rep(1:10, 20)), x = rnorm(n))
d$yb <- pmin(pmax(plogis(d$x / 3 + rnorm(n, 0, .5)), .01), .99)
fit <- suppressWarnings(frm(bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(), data = d))
blk <- vapply(fit$frame$re_blocks, function(b) paste(b$theta_idx, collapse = ","), "")
cat(mode, "theta per block (mu, phi):", blk, "\n")
show <- function(label, pl) {
  e <- frmtmb:::resolve_prior_input(fit, pl)$entries
  cat(sprintf("%-28s %s\n", label, paste(vapply(e, function(z)
    paste0(z$comp, z$idx, "=", z$dist$kind, "(", paste(unlist(z$dist[-1]), collapse = ","), ")"), ""),
    collapse = "  ")))
  if (mode == "lane") {
    t <- validate_prior(pl, bf(yb ~ x + (1 | g), phi ~ (1 | g)) + Beta(), data = d)
    t <- t[t$class == "sd", c("prior", "group", "dpar", "source")]
    cat("   table:", paste(apply(t, 1, paste, collapse = "/"), collapse = " | "), "\n")
  }
}
show("sd THEN sd dpar phi", set_prior("normal(0, 1)", class = "sd") + set_prior("normal(0, 5)", class = "sd", dpar = "phi"))
show("sd dpar phi THEN sd", set_prior("normal(0, 5)", class = "sd", dpar = "phi") + set_prior("normal(0, 1)", class = "sd"))
show("sd group g only", set_prior("normal(0, 5)", class = "sd", group = "g"))
