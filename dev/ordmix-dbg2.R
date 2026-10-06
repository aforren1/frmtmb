.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261005 + 11)
n <- 400
x <- rnorm(n); z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x, z, g)
ob <- frm(bf(y ~ x), family = mixture(cumulative(), cumulative(), cumulative()), data = d, dry_run = "objective")
print(names(ob))
obj <- ob$obj %||% ob
p <- obj$par
print(p)
print(obj$fn(p)); print(obj$gr(p))
tr <- list()
fn <- function(q) { v <- obj$fn(q); tr[[length(tr) + 1]] <<- q; v }
gr <- function(q) { g <- obj$gr(q); if (any(!is.finite(g))) { cat("NaN grad at\n"); print(q); print(g); print(obj$fn(q)) }; g }
o <- try(nlminb(p, fn, gr))
print(o)
