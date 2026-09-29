# Reviewer, claim 1, part 2:
#  (a) the no-smooth identity across builds, by identical() on the saved
#      fitted() arrays;
#  (b) which ROUTE a gaussian distributional fit's Est.Error takes, asked
#      of is.matrix(fitted_point()) rather than of as.matrix();
#  (c) the wall-clock cost the extra differenced coefficients add.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb from:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

# (b) the route
set.seed(53)
nE <- 300
dE <- data.frame(x = stats::runif(nE, -2, 2))
dE$y <- sin(dE$x) + stats::rnorm(nE, 0, exp(-0.5 + 0.4 * sin(2 * dE$x)))
fE <- suppressWarnings(frm(bf(y ~ s(x, k = 6), sigma ~ s(x, k = 6)),
                           family = gaussian(), data = dE))
ndE <- data.frame(x = c(-1.5, 0, 1.5))
p <- get("fitted_point", envir = ns)(fE, ndE, NA)
cat("E gaussian: class(fitted_point) =", paste(class(p), collapse = "/"),
    "| is.matrix =", is.matrix(p), "| length =", length(p), "\n")
cat("E so fitted_point_se takes the", if (is.matrix(p)) "MATRIX (fd)" else
  "SCALAR (analytic)", "route\n")
p2 <- get("fitted_point", envir = ns)(fE, ndE, NA, "response", NULL, "sigma")
cat("E dpar = sigma: is.matrix =", is.matrix(p2), "\n")
saveRDS(list(mu = fitted(fE, newdata = ndE, re_formula = NA),
             sg = fitted(fE, newdata = ndE, re_formula = NA,
                         dpar = "sigma")),
        file.path("dev", paste0("resmooth-rev-dpar-", ARM, ".rds")))

# the near-linear smooth, saved bitwise
set.seed(47)
nI <- 400
dI <- data.frame(x = stats::runif(nI, -2, 2))
latI <- 0.9 * dI$x + stats::rlogis(nI)
dI$y <- factor(cut(latI, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fI <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                           data = dI))
saveRDS(fitted(fI, newdata = data.frame(x = c(-1.5, 0, 1.5)),
               re_formula = NA),
        file.path("dev", paste0("resmooth-rev-nearlin-", ARM, ".rds")))

# (c) cost: the gp() fit with 160 coefficients, and the plain s(x) one
set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
tm <- function(expr, reps = 3) {
  best <- Inf
  for (i in seq_len(reps)) {
    t0 <- proc.time(); force(expr); e <- (proc.time() - t0)[["elapsed"]]
    best <- min(best, e)
  }
  best
}
cat("D gp(x), n_b =", length(fD$estimates[["b"]]),
    "| fitted(newdata, NA) min of 3, seconds:",
    sprintf("%.2f", tm(fitted(fD, newdata = data.frame(x = c(-1.5, 0, 1.5)),
                              re_formula = NA))), "\n")
set.seed(31)
ng <- 40; per <- 20
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- 1.0 * sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
cat("B2 s(x)+(1|g), 40 levels, n_b =", length(fB$estimates[["b"]]),
    "| fitted() in sample, NULL, min of 3, seconds:",
    sprintf("%.2f", tm(fitted(fB, re_formula = NULL))), "\n")
cat("B2 | fitted() in sample, NA, min of 3, seconds:",
    sprintf("%.2f", tm(fitted(fB, re_formula = NA))), "\n")
cat("DONE\n")
