# Review of lane refits: bernoulli reading of the number 1 by refit() vs
# residuals(newdata =) on a fit coded 1/2; gr(cov =, dist = "student")
# under ||.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
set.seed(9)
n <- 40
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.3 + x))
d <- data.frame(x = x, y = s + 1)      # failure 1, success 2
f12 <- frm(y ~ x, data = d, family = bernoulli())
d01 <- data.frame(x = x, y = s)
f01 <- frm(y ~ x, data = d01, family = bernoulli())
nd <- d; nd$y <- 1
go("residuals(f12, newdata y = 1) == residuals(f01, newdata y = 0)?",
   { a <- residuals(f12, newdata = nd, type = "response")[, 1]
     b0 <- residuals(f01, newdata = transform(d01, y = 0), type = "response")[, 1]
     b1 <- residuals(f01, newdata = transform(d01, y = 1), type = "response")[, 1]
     c(reads_1_as_failure = isTRUE(all.equal(a, b0)),
       reads_1_as_success = isTRUE(all.equal(a, b1))) })
go("refit(f12, rep(1, n)) == refit(f01, rep(1, n))? (1 read as success)",
   suppressWarnings(identical(fixef(refit(f12, rep(1, n))),
                              fixef(refit(f01, rep(1, n))))))
go("refit(f12, rep(1, n)) == refit(f01, rep(0, n))? (1 read as failure)",
   suppressWarnings(identical(fixef(refit(f12, rep(1, n))),
                              fixef(refit(f01, rep(0, n))))))
# one row differs: 1 is then a value
v <- rep(1, n); v[1] <- 2
go("refit(f12, c(2, 1, 1, ...)): the 1s read as failures",
   suppressWarnings(identical(fixef(refit(f12, v)),
                              fixef(refit(f01, v - 1)))))

cat("\n==== gr(cov =, dist = 'student') under || ====\n")
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
dd <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                            levels = paste0("g", 1:ng)),
                 f = factor(rep(c("a", "b", "c"), 40)), x = rnorm(120))
dd$y <- rnorm(120) + rnorm(ng)[as.integer(dd$g)]
go("numeric (1 + x || gr(g, cov = A, dist = 'student'))",
   { f <- frm(y ~ x + (1 + x || gr(g, cov = A, dist = "student")), data = dd,
              data2 = list(A = A))
     c(logLik = as.numeric(logLik(f)), n_blocks = length(f$frame$re_blocks)) })
go("factor (0 + f || gr(g, cov = A, dist = 'student'))",
   { f <- frm(y ~ x + (0 + f || gr(g, cov = A, dist = "student")), data = dd,
              data2 = list(A = A))
     c(logLik = as.numeric(logLik(f)), n_blocks = length(f$frame$re_blocks)) })
go("(1 + x || gr(g, cov = A)) block labels",
   { f <- frm(y ~ x + (1 + x || gr(g, cov = A)), data = dd, data2 = list(A = A))
     vapply(f$frame$re_blocks, function(b) b$label %||% NA_character_, "") })
