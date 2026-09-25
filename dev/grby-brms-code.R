# brms 2.23.0's Stan code, Stan data and default priors for by-split
# terms, and its refusal of a level in two by-levels (make_stancode(),
# make_standata() and get_prior(); nothing is compiled).
# Run: Rscript dev/grby-brms-code.R > dev/grby-log/brms-code.txt
.libPaths(c("/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(brms))
set.seed(1)
ng <- 12
d <- data.frame(g = factor(rep(1:ng, each = 5)))
d$f <- factor(ifelse(as.integer(d$g) <= 6, "a", "b"))
d$x <- rnorm(nrow(d)); d$y <- rnorm(nrow(d))
d$g1 <- d$g; d$g2 <- factor(sample(1:ng, nrow(d), TRUE), levels = 1:ng)
d$f1 <- d$f; d$f2 <- factor(ifelse(as.integer(d$g2) <= 6, "a", "b"))
A <- diag(ng); dimnames(A) <- list(levels(d$g), levels(d$g))
cat("==== intercept by\n")
cat(make_stancode(y ~ x + (1 | gr(g, by = f)), data = d))
cat("==== slope by\n")
sc <- make_stancode(y ~ x + (1 + x | gr(g, by = f)), data = d)
cat(sc)
cat("==== prior\n")
print(get_prior(y ~ x + (1 + x | gr(g, by = f)), data = d))
cat("==== cov by\n")
msg <- function(e) conditionMessage(e)
sc <- tryCatch(make_stancode(y ~ x + (1 + x | gr(g, by = f, cov = A)),
                             data = d, data2 = list(A = A)), error = msg)
cat(sc)
cat("==== mm by\n")
sc <- tryCatch(make_stancode(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))),
                             data = d), error = msg)
cat(sc)
print(tryCatch(get_prior(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2))),
                         data = d), error = msg))
cat("==== bad\n")
d2 <- d; d2$f[1] <- "b"
print(tryCatch(make_standata(y ~ x + (1 | gr(g, by = f)), data = d2),
               error = msg))
cat("==== standata\n")
sdat <- make_standata(y ~ x + (1 + x | gr(g, by = f)), data = d)
str(sdat[grep("_1$", names(sdat))])
