# Punch round 2, m1 follow-up: on mo() seed 12, how far does each
# conditional_effects("income:age") row move along the lost zeta2_2
# direction? The cosine jc_nonest() reads, per grid row.
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(12)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
d <- data.frame(income, ls)
d$age <- rnorm(100, mean = 40, sd = 10)
fit <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
jc <- ns$get_joint_cov(fit)
cat("null basis columns:", NCOL(jc$null), "\n")
print(round(jc$null, 4))
nd <- expand.grid(income = factor(lev, levels = lev, ordered = TRUE),
                  age = c(30, 40, 50))
nd$income <- factor(as.character(nd$income), levels = lev, ordered = TRUE)
lb <- frm_lp_basis(fit, newdata = nd)
A <- as.matrix(lb$A)
u <- jc$units %||% rep(1, NROW(jc$null))
Nq <- jc$null / u
Gq <- A * rep(u[lb$coef_pos], each = nrow(A))
cosv <- abs(Gq %*% Nq[lb$coef_pos, , drop = FALSE]) /
  outer(sqrt(rowSums(Gq^2)), sqrt(colSums(Nq^2)))
print(cbind(nd, cos = signif(cosv[, 1], 3), nonest = lb$se_nonest))
cat("coef_pos names:", jc$names[lb$coef_pos], "\n")
