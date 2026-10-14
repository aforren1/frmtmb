# Review round 2: gev_fit_check() (warns when a fitted xi is within 1e-3
# of -1). Seeds 621 and 622 again, and the false-alarm rate on legal
# models with true xi -0.9 to -0.99: a warning on a fit whose maximum is
# interior is a false alarm. Arg 1: link ("log1p" or "identity").
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
link <- commandArgs(TRUE)[1]
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
fitone <- function(y) {
  ws <- character(0)
  f <- withCallingHandlers(
    frm(y ~ 1, family = gen_extreme_value(link_xi = link),
        data = data.frame(y = y)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  xi <- frm_linpred(f, dpar = "xi", type = "response")[1]
  # d logLik / d xi on the natural scale, from the objective's gradient
  gr <- f$obj$gr(f$opt$par)
  ixi <- which(names(f$opt$par) == "betad")[2L]
  dxi <- if (link == "log1p") 1 + xi else 1
  list(xi = xi, dll = -gr[ixi] / dxi, conv = f$opt$convergence,
       edge_warn = any(grepl("the edge at -1", ws)),
       other = any(!grepl("the edge at -1", ws)), ll = as.numeric(logLik(f)))
}
cat("link", link, "\n")
if (link == "identity") {
  for (s in c(621, 622)) {
    set.seed(s)
    y <- rg(60, 0, 1, if (s == 621) -0.8 else -0.95)
    r <- fitone(y)
    cat(sprintf("seed %d: xi + 1 = %.3e conv %d logLik %.4f edge warning %s\n",
                s, r$xi + 1, r$conv, r$ll, r$edge_warn))
  }
}
res <- NULL
for (xi0 in c(-0.9, -0.95, -0.99)) for (n in c(50, 200, 1000)) {
  for (s in 1:40) {
    set.seed(30000 + 1000 * n + s + round(-xi0 * 100))
    r <- fitone(rg(n, 0, 1, xi0))
    res <- rbind(res, data.frame(xi0 = xi0, n = n, seed = s, xi = r$xi,
                                 dll = r$dll, conv = r$conv,
                                 warn = r$edge_warn, other = r$other))
  }
}
# interior: xi above -1 + 1e-6 and a flat log-likelihood in xi there
res$interior <- res$xi > -1 + 1e-6 & abs(res$dll) < 1e-2
agg <- aggregate(cbind(fits = 1, edge = !interior, warned = warn,
                       false_alarm = warn & interior,
                       missed = !warn & !interior,
                       band = xi > -1 + 1e-6 & xi < -0.999) ~ xi0 + n,
                 data = res, FUN = sum)
print(agg, row.names = FALSE)
cat("\nfits with xi in (-1 + 1e-6, -0.999):\n")
print(res[res$xi > -1 + 1e-6 & res$xi < -0.999, ], row.names = FALSE)
cat("\nedge fits, xi + 1 range:", format(range(res$xi[!res$interior] + 1)), "\n")
cat("interior fits, smallest xi:", min(res$xi[res$interior]), "\n")
