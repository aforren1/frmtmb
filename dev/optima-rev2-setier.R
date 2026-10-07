# Reviewer of lane optima, re-check (overlap): test-optima.R's "a lost
# direction's small loading stays out of a kept row" against lane
# setier's current se_tier3() (its R/se-check.R sourced read-only into
# a child of the lane's namespace), to see whether the two lanes'
# changes compose.
#   Rscript dev/optima-rev2-setier.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
st <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-setier/R/se-check.R",
           envir = st)
set.seed(1)
d <- data.frame(x = rnorm(200), z = rnorm(200))
d$y <- 1 + 0.5 * d$x + 0.3 * d$z + rnorm(200)
fit <- frm(y ~ x + z, data = d)
nm <- frmtmb:::outer_par_names(fit)
p <- fit$opt$par
D <- sqrt(diag(fit$obj$he(p)))
np <- length(p)
ix <- match("x", sub("^.*_", "", nm))
iz <- match("z", sub("^.*_", "", nm))
v1 <- numeric(np)
v1[iz] <- 1
v1[ix] <- 0.022
v1 <- v1 / sqrt(sum(v1^2))
Q <- qr.Q(qr(cbind(v1, diag(np))))[, seq_len(np)]
if (sum(Q[, 1] * v1) < 0) Q[, 1] <- -Q[, 1]
S <- Q %*% diag(c(-0.5, rep(1, np - 1))) %*% t(Q)
for (who in c("optima", "setier")) {
  f3 <- if (who == "optima") frmtmb:::se_tier3 else st$se_tier3
  r <- tryCatch(f3(fit, S * outer(D, D), matrix(0, np, np), p,
                   exact = TRUE), error = function(e) e)
  if (inherits(r, "error")) {
    cat(who, "ERROR", conditionMessage(r), "\n")
    next
  }
  cat(sprintf("%-6s lost %s | z lost %s, x kept %s, null row x max %g\n",
              who, paste(names(r$lost), r$lost, sep = "=", collapse = ","),
              nm[iz] %in% names(r$lost), !(nm[ix] %in% names(r$lost)),
              max(abs(r$null[ix, ]))))
}
