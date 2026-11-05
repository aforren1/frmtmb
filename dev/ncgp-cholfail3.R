# Walks frm()'s tape from the factor outward to find where a failed factor stops
# being NaN. Seed 20261102.
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
n <- 60L
set.seed(20261102)
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
y <- sin(2 * x) + 0.1 * e0
d <- data.frame(x = x, y = y)
f2 <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
ns <- asNamespace("frmtmb")
unlockBinding("gp_nugget", ns); assign("gp_nugget", 1e-17, envir = ns)
fr <- f2$frame
bk <- fr$re_blocks[[1]]
th0 <- f2$estimates$theta
thb <- c(th0[1], log(5))
z <- f2$estimates$b
show <- function(tag, f) {
  r <- tryCatch(f(), error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("NCGPC3 %-28s %s\n", tag, paste(format(r, digits = 7), collapse = " ")))
}
t1 <- MakeTape(function(th) frmtmb:::gp_factor(th, bk)[60, 60], th0)
show("gp_factor corner", function() c(t1(th0), t1(thb)))
t2 <- MakeTape(function(th) sum(frmtmb:::gp_field(z, th, bk)), th0)
show("gp_field sum", function() c(t2(th0), t2(thb)))
t3 <- MakeTape(function(th) sum(frmtmb:::expand_b(fr, z, th)), th0)
show("expand_b sum", function() c(t3(th0), t3(thb)))
nll <- frmtmb:::build_objective(fr)
p0 <- f2$estimates
t4 <- MakeTape(function(a) { p <- p0; p$theta <- a; nll(p) }, th0)
show("objective, theta only", function() c(t4(th0), t4(thb)))
o <- RTMB::MakeADFun(nll, p0, random = "b", silent = TRUE)
pb <- f2$opt$par; pb[which(names(pb) == "theta")[2]] <- log(5)
show("MakeADFun plain", function() c(o$fn(f2$opt$par), o$fn(pb)))
show("  env f at b = 0", function() { fp <- o$env$last.par; fp[names(fp) == "b"] <- 0; fp[names(fp) == "theta"] <- thb; o$env$f(fp) })
show("  env hess", function() { fp <- o$env$last.par; fp[names(fp) == "theta"] <- thb; h <- o$env$spHess(fp, random = TRUE); range(as.matrix(h)) })
print(body(o$fn))
print(o$env$inner.control)
print(o$env$random.start)
