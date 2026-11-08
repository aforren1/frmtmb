# Review of lane readers: print the pooled lognormal and poisson bands of
# both arms at five grid points (from dev/review-readers-pool.R's RDS).
a <- readRDS("dev/review-readers-out/pool-before.rds")
b <- readRDS("dev/review-readers-out/pool-after.rds")
for (nm in c("poisson", "lognormal", "hurdle")) {
  i <- c(1, 25, 50, 75, 100)
  cat("==", nm, "\n")
  print(data.frame(x = a[[nm]]$x[i], est = a[[nm]]$estimate__[i],
                   lo_before = a[[nm]]$lower__[i], lo_after = b[[nm]]$lower__[i],
                   up_before = a[[nm]]$upper__[i], up_after = b[[nm]]$upper__[i],
                   se_before = a[[nm]]$se__[i], se_after = b[[nm]]$se__[i]),
        digits = 4)
}
