# Reviewer, fix round 1: row_lpdf_prob() (families without cdf_log)
# against 0.72.0. The objective at the start of each model on the same
# data, printed with its start vector; run with ARM=base and ARM=lane
# and compare the two logs (cdfs-rev-legacy-cmp.R).
arm <- Sys.getenv("ARM", "lane")
lp <- c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(31)
n <- 120
x <- rnorm(n)
cc <- sample(c(-1, 0, 0, 1, 2), n, TRUE)
pos <- data.frame(x = x, y = rweibull(n, 1.5, exp(0.3 * x)) + 0.05, cc = cc)
pos$y2 <- pos$y * 1.6
pos$lb <- 0.04; pos$ub <- max(pos$y2) + 1
gev <- data.frame(x = x, y = rnorm(n, 0.3 * x), cc = cc)
gev$y2 <- gev$y + 0.7
dw <- data.frame(x = x, y = rgeom(n, 0.3), cc = cc)
dw$y2 <- dw$y + 2
fams <- list(exponential = list(exponential(), pos), weibull = list(weibull(), pos),
             frechet = list(frechet(), pos), gen_extreme_value = list(gen_extreme_value(), gev),
             discrete_weibull = list(discrete_weibull(), dw))
forms <- list(cens = "y | cens(cc, y2) ~ x", trunc_both = "y | trunc(lb = LB, ub = UB) ~ x",
              trunc_lb = "y | trunc(lb = LB) ~ x", trunc_ub = "y | trunc(ub = UB) ~ x",
              cens_trunc = "y | cens(cc, y2) + trunc(lb = LB, ub = UB) ~ x")
for (nm in names(fams)) {
  fm <- fams[[nm]][[1]]; d <- fams[[nm]][[2]]
  disc <- nm == "discrete_weibull"
  LB <- if (disc) 1 else if (nm == "gen_extreme_value") min(d$y) - 0.5 else 0.04
  UB <- if (disc) max(d$y2) + 3 else max(d$y2) + 1
  if (disc) d <- d[d$y >= 1, ]
  for (fk in names(forms)) {
    f <- as.formula(gsub("UB", format(UB, digits = 17), gsub("LB", format(LB, digits = 17), forms[[fk]])))
    r <- tryCatch({
      ob <- frm(f, family = fm, data = d, dry_run = "objective")$obj
      sprintf("fn %.15g | par %s", ob$fn(ob$par), paste(format(ob$par, digits = 10), collapse = " "))
    }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 100)))
    cat(sprintf("%-18s %-11s %s\n", nm, fk, r))
  }
}
if (arm == "lane") {
  cat("== discrete_weibull, four codes: lane objective vs brms's convention by hand\n")
  d <- dw
  ob <- frm(y | cens(cc, y2) ~ x, family = discrete_weibull(), data = d, dry_run = "objective")$obj
  fam <- frmtmb:::as_frmtmb_family(discrete_weibull())
  p <- ob$par
  # dpars at the start: mu on its logit link, shape on log, by name
  lin <- p[1] + p[2] * d$x
  dp <- list(mu = 1 / (1 + exp(-lin)), shape = exp(p[3]))
  ll <- fam$lpdf(d$y, dp, list())
  Fy <- fam$lcdf(d$y, dp, list()); F2 <- fam$lcdf(d$y2, dp, list())
  ll[d$cc == -1] <- log(Fy[d$cc == -1])
  ll[d$cc == 1] <- log1p(-Fy[d$cc == 1])
  ll[d$cc == 2] <- log(F2[d$cc == 2] - Fy[d$cc == 2])
  cat(sprintf("objective %.12f  hand (brms reading) %.12f  par %s\n", ob$fn(p), -sum(ll),
              paste(signif(p, 6), collapse = " ")))
  print(names(fam$links)); print(sapply(fam$links, function(l) l$name))
}
