# brms_car_const() of the spatial suite on the data of dev/latres-spatial.R, against the L1 constants of its log.
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib","C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spatial)})
hel <- new.env(parent = asNamespace("frmtmb.spatial")); hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)
sys.source("extensions/frmtmb.spatial/tests/testthat/helper-brms.R", envir = hel)
set.seed(167); nl <- 12L; M <- matrix(0, nl, nl)
for (i in seq_len(nl - 1L)) M[i, i + 1L] <- M[i + 1L, i] <- 1
M[1, nl] <- M[nl, 1] <- 1; locs <- paste0("L", seq_len(nl)); dimnames(M) <- list(locs, locs)
d <- data.frame(g = factor(rep(locs, each = 8L), levels = locs)); d$x <- rnorm(nrow(d)); d$count <- 1L
for (ty in c("escar", "icar", "esicar")) {
  f <- as.formula(sprintf("count ~ x + car(M, gr = g, type = '%s')", ty))
  sd <- hel$brms_standata(brms::bf(f), data = d, family = brms::brmsfamily("poisson"), data2 = list(M = M))
  cat(ty, format(hel$brms_car_const(sd, ty), digits = 8), "\n")
}
