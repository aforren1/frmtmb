# Item 3.5, collapse arm: per replicate, lambda's logit, its standard
# error, the log-likelihood gap to the plain family, and whether
# diagnose() named lambda at the end of its link.
x <- lapply(list.files("dev/phase3b-log/recov", "^collapse-[0-9]+[.]rds$",
                       full.names = TRUE), readRDS)
for (i in seq_along(x)) {
  if (!is.null(x[[i]]$error)) {
    f2 <- sprintf("dev/phase3b-log/recov-rerun/collapse-%d.rds", x[[i]]$seed)
    if (file.exists(f2)) x[[i]] <- readRDS(f2)
  }
}
x <- x[vapply(x, function(r) is.null(r$error), TRUE)]
tab <- do.call(rbind, lapply(x, function(r) {
  ci <- r$confint["lambda_(Intercept)", ]
  data.frame(seed = r$seed, eta = ci[["est"]],
             se = (ci[["upr"]] - ci[["lwr"]]) / 3.92,
             dll = r$loglik - (if (is.null(r$loglik_plain)) NA else r$loglik_plain),
             named = any(grepl("end of its link: lambda", r$diagnose)))
}))
tab <- tab[order(tab$named, tab$eta), ]
print(tab, digits = 4, row.names = FALSE)
cat("not named:", sum(!tab$named), "; their eta range",
    format(range(tab$eta[!tab$named]), digits = 4), "; se range",
    format(range(tab$se[!tab$named]), digits = 4), "\n")
cat("named: eta range", format(range(tab$eta[tab$named]), digits = 4),
    "; se range", format(range(tab$se[tab$named]), digits = 4), "\n")
cat("|logLik gap| max over all:", format(max(abs(tab$dll), na.rm = TRUE)),
    "; NA gaps (plain fit failed):", sum(is.na(tab$dll)), "\n")
