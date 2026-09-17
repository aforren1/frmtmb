## Reviewer recheck: brms 2.23.0 on the collision classes of
## dev/brmsnames-rev2-collide.R. C4 is fitted (fresh compile, chains 1,
## iter 60, seed 1); C2, C3 and C5 are parsed only.
##   Rscript dev/brmsnames-rev2-collide-brms.R
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(packageVersion("StanHeaders") == "2.32.10")
suppressMessages(library(brms))
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
d$fab3 <- factor(ifelse(seq_len(nrow(d)) %% 3 == 0, "0",
                        as.character(d$fab)), levels = c("0", "a b", "ab"))
gd <- as.character(d$g)
i1 <- which(gd == "lvl 1")
gd[i1[c(TRUE, FALSE)]] <- "lvl.1"
d$gd <- factor(gd)
cat("== C2 ==\n")
print(try1(make_standata(y ~ fab3, data = d)$X))
r <- try1(brm(y ~ fab3, data = d, chains = 0))
print(if (is.character(r)) r else "C2: brm(chains = 0) built")
cat("== C3 ==\n")
r <- try1(brm(y ~ I(x^2) + IxE2, data = d, chains = 0))
print(if (is.character(r)) r else "C3: brm(chains = 0) built")
cat("== C5 ==\n")
r <- try1(brm(mvbf(bf(y_a ~ x), bf(ya ~ x), rescor = FALSE), data = d,
              chains = 0))
print(if (is.character(r)) r else "C5: brm(chains = 0) built")
cat("== C4 ==\n")
b <- try1(brm(y ~ x + (1 | gd), data = d, chains = 1, iter = 60,
              warmup = 30, seed = 1, refresh = 0, backend = "rstan"))
if (is.character(b)) print(b) else {
  saveRDS(b, "dev/stan-cache/brmsnames-rev2-brms-C4.rds")
  v <- variables(b)
  print(grep("^r_", v, value = TRUE))
  print(try1(dimnames(ranef(b)$gd)[[1]]))
  print(try1(dim(as_draws_df(b))))
}
