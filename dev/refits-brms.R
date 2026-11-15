# brms 2.23.0 on the lane's parity questions; no Stan compile needed.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  r <- tryCatch(expr, error = function(e)
    paste("ERROR:", gsub("[\r\n]+", " ", conditionMessage(e))))
  print(r)
  invisible(r)
}

## A. character ordinal response
set.seed(2501)
n <- 60
d <- data.frame(x = rnorm(n))
lab <- c("none", "mild", "moderate", "severe")
idx <- pmin(4, pmax(1, round(2.5 + d$x + rnorm(n, 0, .7))))
d$y <- lab[idx]
d$yn <- as.character(idx)
go("A char labels cumulative", {
  sd <- standata(y ~ x, family = cumulative(), data = d); table(sd$Y) })
go("A char digits cumulative", {
  sd <- standata(yn ~ x, family = cumulative(), data = d); table(sd$Y) })
go("A char labels sratio", {
  sd <- standata(y ~ x, family = sratio(), data = d); table(sd$Y) })
d$yu <- factor(d$y)
go("A unordered factor", {
  sd <- standata(yu ~ x, family = cumulative(), data = d); table(sd$Y) })

## C. y ~ . - w: the formula brms stores, and standata on data without w
set.seed(3)
d4 <- data.frame(x = rnorm(80), z = rnorm(80), w = rnorm(80))
d4$y <- d4$x + rnorm(80)
go("C stored formula", {
  f <- brm(y ~ . - w, data = d4, empty = TRUE)
  print(f$formula)
  d5 <- d4; d5$w <- NULL
  u <- tryCatch(update(f, newdata = d5, chains = 0),
                error = function(e) paste("update ERROR:",
                                          conditionMessage(e)))
  print(u)
  standata(f$formula, data = d5)$K
})

## E. (0 + f || gr(g, cov = A)): the Stan model for a factor under ||
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d8 <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                            levels = paste0("g", 1:ng)),
                 f = factor(rep(c("a", "b"), 60)))
d8$y <- rnorm(120)
go("E stancode", {
  sc <- stancode(y ~ 1 + (0 + f || gr(g, cov = A)), data = d8,
                 data2 = list(A = A))
  cat(grep("r_1|sd_1|Lcov|z_1|M_1", strsplit(sc, "\n")[[1]],
           value = TRUE), sep = "\n")
  invisible(NULL)
})
go("E prior table", {
  print(default_prior(y ~ 1 + (0 + f || gr(g, cov = A)), data = d8,
                      data2 = list(A = A)))
  invisible(NULL)
})

## F. (cs(1) | g) with sratio
set.seed(12)
d9 <- data.frame(g = factor(rep(1:15, each = 20)), x = rnorm(300))
d9$y <- sample(1:4, 300, TRUE)
go("F stancode", {
  sc <- stancode(y ~ x + (cs(1) | g), family = sratio(), data = d9)
  cat(grep("r_1|sd_1|L_1|z_1|mucs|J_1|Z_1|M_1|NC_1|lpdf",
           strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
  invisible(NULL)
})
go("F prior table", {
  print(default_prior(y ~ x + (cs(1) | g), family = sratio(), data = d9))
  invisible(NULL)
})
go("F cumulative refused?", {
  stancode(y ~ x + (cs(1) | g), family = cumulative(), data = d9)
  "no error"
})
go("F (cs(1) || g)", {
  sc <- stancode(y ~ x + (cs(1) || g), family = acat(), data = d9)
  cat(grep("r_1|sd_1|L_1|mucs", strsplit(sc, "\n")[[1]], value = TRUE),
      sep = "\n")
  invisible(NULL)
})
go("F (1 + cs(x) | g)", {
  sc <- stancode(y ~ x + (1 + cs(x) | g), family = sratio(), data = d9)
  cat(grep("r_1|sd_1|L_1|mucs|  mu", strsplit(sc, "\n")[[1]], value = TRUE),
      sep = "\n")
  invisible(NULL)
})
cat("\nDONE\n")
