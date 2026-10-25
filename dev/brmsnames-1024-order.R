## Row order of default_prior() and validate_prior(), brms 2.23.0 against
## frmtmb, on the rows both have (frmtmb's own class "theta" left out).
##
##   Rscript dev/brmsnames-1024-order.R base|lane > dev/brmsnames-1024-log/order-<arm>.txt
##
## Data: dev/brmsnames-1024-cases.R (seed 1024). A row is
## class/coef/group/resp/dpar/nlpar; brms's rows that frmtmb refuses by
## class name (sds, meanme, sdme) and frmtmb's that brms lacks are left
## out of both before the order is compared.
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
source("dev/brmsnames-1024-cases.R")
d <- brmsnames_data()
forms <- list(
  list("y ~ exp(x) + z + (1 + z | g)", "gaussian"),
  list("y ~ x + f, sigma ~ z", "gaussian"),
  list("y ~ x + (1 | g), sigma ~ z + (1 | g)", "gaussian"),
  list("yc ~ log(xp) + (1 | g) + (0 + x | f)", "poisson"),
  list("y ~ a + b * x, a ~ exp(z) + (1 | g), b ~ 1 + z, nl = TRUE",
       "gaussian"))
key <- function(p) {
  p <- as.data.frame(p)
  p <- p[!p$class %in% c("theta", "sds", "meanme", "sdme"), ]
  paste(p$class, p$coef, p$group, p$resp, p$dpar, p$nlpar, sep = "/")
}
for (f in forms) {
  bb <- function(pkg) {
    env <- list2env(list(bf = getExportedValue(pkg, "bf")))
    eval(str2lang(paste0("bf(", f[[1]], ")")), env)
  }
  bd <- key(q(brms::default_prior(bb("brms"), data = d, family = f[[2]])))
  fd <- key(q(default_prior(bb("frmtmb"), data = d, family = f[[2]])))
  pr <- "normal(0, 5)"
  nlp <- if (grepl("nl = TRUE", f[[1]], fixed = TRUE)) "a" else ""
  bv <- key(q(brms::validate_prior(brms::set_prior(pr, class = "b",
                                                   nlpar = nlp),
                                   bb("brms"), data = d, family = f[[2]])))
  fv <- key(q(validate_prior(set_prior(pr, class = "b", nlpar = nlp),
                             bb("frmtmb"), data = d, family = f[[2]])))
  fd <- fd[fd %in% bd]
  bd <- bd[bd %in% fd]
  fv <- fv[fv %in% bv]
  bv <- bv[bv %in% fv]
  cat(sprintf("%-60s default_prior same order: %s; validate_prior: %s\n",
              f[[1]], identical(fd, bd), identical(fv, bv)))
  if (!identical(fd, bd)) {
    cat("  brms  :", bd, "\n  frmtmb:", fd, "\n")
  }
  if (!identical(fv, bv)) {
    cat("  brms v  :", bv, "\n  frmtmb v:", fv, "\n")
  }
}
