# Lane postfix: the objective of brms_nonlinear's fit_loss (no priors)
# over the intercept grid, to see where the scan lands and why.
#   Rscript dev/postfix-nlscan-loss.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
data("loss", package = "brms")
f7 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
         ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE)
u <- frm(f7, data = loss, family = gaussian(), dry_run = "objective")
obj <- u$obj
cat("outer:", names(obj$par), "\n")
cat("start:", obj$par, "\n")
f <- function(p) tryCatch(as.numeric(obj$fn(p)), error = function(e) NA)
p0 <- obj$par
pos <- which(names(p0) == "beta")
grid <- c(0, as.vector(t(c(1, -1) %o% 10^(-2:3))))
cat("hand (5000, 1, 45):", f(replace(p0, pos, c(5000, 1, 45))), "\n")
res <- expand.grid(ult = grid, omega = grid, theta = grid)
res$nll <- apply(res, 1, function(r) f(replace(p0, pos, r[1:3])))
res <- res[order(res$nll), ]
print(head(res, 12))
g125 <- c(0, as.vector(t(c(1, -1) %o% as.vector(c(1, 2, 5) %o%
                                                  10^(-2:3)))))
res2 <- expand.grid(ult = g125, omega = g125, theta = g125)
res2$nll <- apply(res2, 1, function(r) f(replace(p0, pos, r[1:3])))
res2 <- res2[order(res2$nll), ]
cat("1-2-5 grid, best 8:\n")
print(head(res2, 8))
# what nlminb reaches from the best 1-2-5 point and from the best
# power-of-ten point, holding nothing
for (st in list(unlist(res[1, 1:3]), unlist(res2[1, 1:3]))) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    frm(f7, data = loss, family = gaussian(),
        start = list(beta = unname(st))))),
    error = function(e) conditionMessage(e))
  cat("from", st, ":", if (is.character(fit)) fit else
    format(as.numeric(logLik(fit)), digits = 10), "\n")
}
cat("common values:\n")
for (cv in grid) cat(sprintf("  c = %8g  nll = %s\n", cv,
                             format(f(replace(p0, pos, cv)))))
