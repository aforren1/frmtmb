# Review of lane refits: brms 2.23.0 behavior (standata, no compile).
# 2939 group counts, 1259 character ordinal, 1567 dot formulas,
# bernoulli newdata coding.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  r <- tryCatch(expr, error = function(e)
    paste0("ERROR: ", conditionMessage(e)))
  print(r)
  invisible(r)
}
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                           levels = paste0("g", 1:ng)),
                f = factor(rep(c("a", "b", "c"), 40)),
                f2 = factor(rep(c("a", "b"), 60)),
                f4 = factor(rep(c("a", "b", "c", "d"), 30)),
                h = factor(rep(c("u", "v"), each = 60)),
                x = rnorm(120), z = rnorm(120), y = rnorm(120))
forms <- list(
  num_cov = y ~ x + (1 + x || gr(g, cov = A)),
  f_cov = y ~ x + (1 + f || gr(g, cov = A)),
  f0_cov = y ~ x + (0 + f || gr(g, cov = A)),
  f2_cov = y ~ x + (1 + f2 || gr(g, cov = A)),
  f4_0_cov = y ~ x + (0 + f4 || gr(g, cov = A)),
  mix_cov = y ~ x + (1 + z + f || gr(g, cov = A)),
  by_num_cov = y ~ x + (1 + x || gr(g, cov = A, by = h)),
  by_f_cov = y ~ x + (0 + f2 || gr(g, cov = A, by = h)),
  plain_f = y ~ x + (1 + f || g)
)
for (nm in names(forms)) {
  go(nm, {
    sd <- standata(forms[[nm]], data = d, data2 = list(A = A))
    dp <- default_prior(forms[[nm]], data = d, data2 = list(A = A))
    list(M = unlist(sd[grep("^M_", names(sd))]),
         N = unlist(sd[grep("^N_[0-9]", names(sd))]),
         NC = unlist(sd[grep("^NC_", names(sd))]),
         Z = grep("^Z_", names(sd), value = TRUE),
         sd_rows = paste(dp$group, dp$coef)[dp$class == "sd"],
         Jby = unlist(sd[grep("^Jby_|^Nby_", names(sd))]))
  })
}
go("stancode by_num_cov r_1 lines", {
  sc <- stancode(forms$by_num_cov, data = d, data2 = list(A = A))
  grep("r_1|sd_1", strsplit(sc, "\n")[[1]], value = TRUE)
})

cat("\n==== 1259 character ordinal ====\n")
d2 <- data.frame(x = rnorm(30))
d2$yc <- sample(c("none", "mild", "severe"), 30, TRUE)
d2$yd <- sample(c("1", "2", "3"), 30, TRUE)
d2$yn <- as.numeric(d2$yd)
go("cumulative character labels",
   standata(yc ~ x, data = d2, family = cumulative())$Y)
go("cumulative character digits",
   standata(yd ~ x, data = d2, family = cumulative())$Y)
go("cumulative numeric", head(standata(yn ~ x, data = d2,
                                      family = cumulative())$Y))
go("sratio character digits",
   standata(yd ~ x, data = d2, family = sratio())$Y)

cat("\n==== 1567 dot formulas ====\n")
d3 <- data.frame(y = rnorm(20), x = rnorm(20), z = rnorm(20), w = rnorm(20),
                 o = abs(rnorm(20)))
d3n <- d3[, c("y", "x", "z", "o")]
go("y ~ . - w, with w", colnames(standata(y ~ . - w, data = d3)$X))
go("y ~ . - w, without w", colnames(standata(y ~ . - w, data = d3n)$X))
go("y ~ .^2 - w, with w", colnames(standata(y ~ .^2 - w, data = d3)$X))
go("y ~ .^2 - w, without w", colnames(standata(y ~ .^2 - w, data = d3n)$X))
go("y ~ . - w - o + offset(o), with w",
   { s <- standata(y ~ . - w - o + offset(o), data = d3)
     list(colnames(s$X), head(s$offsets)) })
go("y ~ . - w - o + offset(o), without w",
   { s <- standata(y ~ . - w - o + offset(o), data = d3n)
     list(colnames(s$X), head(s$offsets)) })
go("y ~ . - q (q absent)", colnames(standata(y ~ . - q, data = d3)$X))

cat("\n==== bernoulli coding of a response ====\n")
d4 <- data.frame(x = rnorm(10))
for (yy in list(c(1, 2, 1, 2, 2, 1, 1, 2, 2, 2), rep(1, 10), rep(2, 10),
                rep(0, 10), c(-1, 1, -1, 1, 1, -1, 1, 1, -1, 1))) {
  d4$y <- yy
  go(paste("bernoulli y =", paste(unique(yy), collapse = ",")),
     standata(y ~ x, data = d4, family = bernoulli())$Y)
}
