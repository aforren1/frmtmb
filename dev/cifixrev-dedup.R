# Reviewer: is "one position" decided exactly by pred_design()'s and
# gp_krig_cov()'s dedup? Usage: Rscript dev/cifixrev-dedup.R <lib|base>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(17)
xg <- round(runif(120, 0, 10), 1)
dg <- data.frame(y = sin(xg) + rnorm(120, 0, 0.3), x = xg,
                 f = factor(sample(c("a", "b"), 120, TRUE)),
                 x2 = runif(120, 0, 5))
fg <- suppressWarnings(frm(bf(y ~ gp(x)), data = dg))
x1 <- 1 / 3
x2 <- x1 * (1 + 2^-52)
x3 <- 1 / 3 + 1e-14
cat("x1 == x2:", x1 == x2, " as.character equal:",
    as.character(x1) == as.character(x2), "\n")
nd <- data.frame(x = c(x1, x2, x3, 5.55, 5.55))
lb <- frm_lp_basis(fg, newdata = nd, extra_cov = TRUE)
A <- unname(as.matrix(lb$A))
E <- unname(as.matrix(lb$extra_cov))
cat("rows 1,2 (distinct doubles 1 ulp apart) A identical:",
    identical(A[1, ], A[2, ]), " max |dA|:", max(abs(A[1, ] - A[2, ])),
    "\n")
cat("rows 1,3 (1e-14 apart) A identical:", identical(A[1, ], A[3, ]),
    " max |dA|:", max(abs(A[1, ] - A[3, ])), "\n")
cat("E[1,2] == E[1,1]:", E[1, 2] == E[1, 1], " E[1,1] - E[1,2]:",
    E[1, 1] - E[1, 2], "\n")
cat("rows 4,5 (one position) A identical:", identical(A[4, ], A[5, ]),
    " E[4,5] == E[4,4]:", E[4, 5] == E[4, 4], "\n")
# gp(x, by = f): the same x in two levels is two sub-GPs
fb <- suppressWarnings(frm(bf(y ~ f + gp(x, by = f)), data = dg))
ndb <- data.frame(x = c(5.55, 5.55, 5.55), f = factor(c("a", "b", "a"),
                                                      levels = c("a", "b")))
lbb <- frm_lp_basis(fb, newdata = ndb, extra_cov = TRUE)
Ab <- unname(as.matrix(lbb$A))
Eb <- unname(as.matrix(lbb$extra_cov))
cat("by = f: rows a and b load disjoint columns:",
    all(Ab[1, ] * Ab[2, ] == 0), " cov(a, b):", Eb[1, 2],
    " rows a and a identical:", identical(Ab[1, ], Ab[3, ]),
    " cov(a, a) == var(a):", Eb[1, 3] == Eb[1, 1], "\n")
# a by level the fit never saw
ndn <- data.frame(x = 5.55, f = factor("c", levels = c("a", "b", "c")))
r <- tryCatch(frm_lp_basis(fb, newdata = ndn, allow_new_levels = TRUE),
              error = function(e) conditionMessage(e))
cat("new by level:", if (is.character(r)) substr(r, 1, 150) else
  paste("extra_var", r$extra_var), "\n")
# two-dimensional positions that agree in one coordinate only
f2 <- suppressWarnings(frm(bf(y ~ gp(x, x2)), data = dg))
nd2 <- data.frame(x = c(5.55, 5.55, 5.55), x2 = c(1.01, 1.02, 1.01))
l2 <- frm_lp_basis(f2, newdata = nd2, extra_cov = TRUE)
A2 <- unname(as.matrix(l2$A))
E2 <- unname(as.matrix(l2$extra_cov))
cat("2-D: rows 1,2 differ:", !identical(A2[1, ], A2[2, ]),
    " rows 1,3 identical:", identical(A2[1, ], A2[3, ]),
    " E[1,2] < E[1,1]:", E2[1, 2] < E2[1, 1], "\n")
# cost: a 2000-row grid with every position repeated once
ndl <- data.frame(x = rep(seq(0.01, 9.99, length.out = 1000), 2))
t1 <- system.time(for (i in 1:3) frm_lp_basis(fg, newdata = ndl,
                                               extra_cov = TRUE))
ndu <- data.frame(x = seq(0.01, 9.99, length.out = 2000))
t2 <- system.time(for (i in 1:3) frm_lp_basis(fg, newdata = ndu,
                                               extra_cov = TRUE))
cat("3 x frm_lp_basis(extra_cov) 2000 rows, 1000 positions (s):",
    t1[["elapsed"]], "; 2000 distinct (s):", t2[["elapsed"]], "\n")
