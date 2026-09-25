# Item 3.5: do the `cont` replicates that stopped at false convergence
# differ from the ones that did not? Estimates and interval hits by
# convergence code.
x <- lapply(list.files("dev/phase3b-log/recov", "^cont-[0-9]+[.]rds$",
                       full.names = TRUE), readRDS)
# an errored first-pass replicate is replaced by its rerun, as in
# dev/phase3b-summarise.R
for (i in seq_along(x)) {
  if (!is.null(x[[i]]$error)) {
    f2 <- sprintf("dev/phase3b-log/recov-rerun/cont-%d.rds", x[[i]]$seed)
    if (file.exists(f2)) x[[i]] <- readRDS(f2)
  }
}
x <- x[vapply(x, function(r) is.null(r$error), TRUE)]
row <- function(r, nm, truth) {
  if (is.null(r$confint) || !nm %in% rownames(r$confint)) return(c(est = NA, hit = NA, se = NA))
  ci <- r$confint[nm, ]
  c(est = ci[["est"]], hit = ci[["lwr"]] < truth && truth < ci[["upr"]],
    se = (ci[["upr"]] - ci[["lwr"]]) / 3.92)
}
tab <- do.call(rbind, lapply(x, function(r) {
  data.frame(seed = r$seed, code = r$conv,
             lam = row(r, "lambda_(Intercept)", qlogis(0.05))[["est"]],
             lam_hit = row(r, "lambda_(Intercept)", qlogis(0.05))[["hit"]],
             lam_se = row(r, "lambda_(Intercept)", qlogis(0.05))[["se"]],
             ndt = row(r, "ndt_(Intercept)", 0)[["est"]],
             cond = row(r, "cond", 0.9)[["est"]],
             cond_hit = row(r, "cond", 0.9)[["hit"]],
             ncont = r$n_cont)
}))
print(tab[order(tab$code, tab$seed), ], digits = 4)
print(aggregate(cbind(lam, lam_hit, lam_se, ndt, cond, cond_hit) ~ code,
                data = tab, FUN = mean), digits = 4)
ok <- is.finite(tab$lam)
cat(sprintf("logit lambda: sd over %d fits %.4f, mean Wald se %.4f, ratio %.2f\n",
            sum(ok), sd(tab$lam[ok]), mean(tab$lam_se[ok]),
            sd(tab$lam[ok]) / mean(tab$lam_se[ok])))
dev <- tab$lam[ok] - qlogis(tab$ncont[ok] / 12000)
cat(sprintf("estimate minus realized share (logit): mean %.4f, sd %.4f\n",
            mean(dev), sd(dev)))
hi <- vapply(x, function(r) r$crange[2], 0)
lo <- vapply(x, function(r) r$crange[1], 0)
cat("observed range top: ", format(summary(hi), digits = 4), "\n")
cat("observed range bottom: ", format(summary(lo), digits = 4), "\n")
cat(sprintf("cor(log width of range, lambda - realized) = %.3f\n",
            cor(log(hi - lo)[ok], dev)))
cat(sprintf("after the width correction log(width / 4.9): mean %.4f, sd %.4f\n",
            mean(dev - log((hi - lo)[ok] / 4.9)), sd(dev - log((hi - lo)[ok] / 4.9))))
