# Lane setier: deg-6 raw polynomial + (1 | g) of dev/setier-collin.R:
# the finite-difference Hessian's spectrum against its own noise.
args <- "src"
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-setier",
                                   quiet = TRUE))
deg <- 6
set.seed(deg)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
fo <- stats::as.formula(paste0("y ~ ", paste0("I(x^", 1:deg, ")",
                                             collapse = " + "), " + (1 | g)"))
f <- suppressWarnings(frm(fo, data = d))
h <- f$cache$hessian_fixed
H <- h$H; E <- h$E
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
Es <- E / outer(D, D)
nk <- sqrt(colSums((Es %*% e$vectors)^2))
print(rbind(ev = e$values / max(e$values), nk = nk / max(e$values)))
print(sdr_of(f)$se_lost)
cat("rcond", rcond(H), "\n")
