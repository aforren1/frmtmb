# Does brms refuse the animal-model design (1 | gr(id, cov = A)) + (1 | id)?
# frmtmb's own suite fits it (test-review-v29.R), and it is identified
# because the two blocks have different covariance matrices over levels.
suppressMessages(library(brms))
set.seed(63)
ng <- 10
A <- diag(ng); A[1, 2] <- A[2, 1] <- 0.5
rownames(A) <- colnames(A) <- paste0("i", 1:ng)
dd <- data.frame(id = factor(rep(rownames(A), each = 4)), y = rnorm(40))
dd$id2 <- dd$id
f <- function(e) tryCatch({force(e); "accepted"}, error = function(err) conditionMessage(err))
cat("brms (1|gr(id, cov = A)) + (1|id):  ",
    f(default_prior(y ~ (1 | gr(id, cov = A)) + (1 | id), dd, data2 = list(A = A))), "\n")
cat("brms (1|gr(id, cov = A)) + (1|id2): ",
    f(default_prior(y ~ (1 | gr(id, cov = A)) + (1 | id2), dd, data2 = list(A = A))), "\n")
