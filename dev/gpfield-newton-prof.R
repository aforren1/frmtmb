# Lane gpfield fix round: where gp_field_newton()'s time goes, Rprof over
# 100 replicates of dev/gpfield-newton.R's poisson n 300 fit (set.seed(5)).
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
ns <- asNamespace("frmtmb")
n <- 300
set.seed(9300 + n)
d <- data.frame(x = sort(round(stats::runif(n, 0, 50), 4)))
d$y <- stats::rpois(n, exp(1 + sin(d$x / 3)))
fit <- gf_fit(bf(y ~ gp(x)), data = d, family = poisson())$fit
map <- ns$outer_par_map(fit)
v0 <- ns$fit_outer_vector(fit, map)
R <- chol(vcov(fit, full = TRUE))
plan <- ns$gp_field_draw_plan(fit, ns$gp_drawn_b(fit))
set.seed(5)
fss <- lapply(1:100, function(s) {
  ns$fit_set_outer(fit, v0 + as.vector(t(R) %*% rnorm(length(v0))), map)
})
steps <- 0L
tr <- ns$gp_field_newton
Rprof(file.path(GF_OUT, "newton-prof.out"), interval = 0.005)
t0 <- proc.time()[[3]]
for (fs in fss) ns$gp_field_newton(fs, plan$newton)
cat("NEWTONPROF 100 replicates", proc.time()[[3]] - t0, "s\n")
Rprof(NULL)
s <- summaryRprof(file.path(GF_OUT, "newton-prof.out"))
print(head(s$by.total, 25))
