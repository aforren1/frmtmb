# Item 9: non-gp fits and gp(k =) digit for digit against rellib-r16. Writes one line per
# reading, in %a, and counts them: the comparison must not pass on two empty files.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
ns <- asNamespace("frmtmb")
set.seed(9009)
n <- 120
d <- data.frame(x = runif(n, 0, 6), z = rnorm(n), g = factor(sample(letters[1:8], n, TRUE)), h = factor(sample(LETTERS[1:5], n, TRUE)),
                o = ordered(sample(1:4, n, TRUE)), t = rep(1:12, 10), id = factor(rep(1:10, each = 12)))
d$y <- 1 + sin(d$x) + 0.5 * d$z + rnorm(8, sd = .6)[as.integer(d$g)] + rnorm(n, sd = 0.4)
d$y2 <- cos(d$x) + rnorm(n, sd = 0.3)
d$cnt <- rpois(n, exp(1 + 0.4 * sin(d$x))); d$bin <- rbinom(n, 1, plogis(sin(d$x))); d$tr <- 10; d$bn <- rbinom(n, 10, plogis(sin(d$x)))
d$pos <- exp(0.3 * d$y) + 0.01; d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$zm <- d$z; d$zm[seq(5, n, 15)] <- NA; d$zobs <- d$z + rnorm(n, sd = 0.1); d$sdz <- 0.1
sp <- list(
  lm = list(y ~ z + x), glmm = list(y ~ z + (1 | g)), slope = list(y ~ z + (z | g) + (1 | h)), smooth = list(y ~ s(x) + z), t2 = list(y ~ t2(x, z)),
  pois = list(cnt ~ z + (1 | g), family = poisson()), bern = list(bin ~ z + (1 | g), family = bernoulli()), binom = list(bn | trials(tr) ~ z, family = binomial()),
  gamma = list(pos ~ z, family = Gamma(link = "log")), student = list(y ~ z, family = student()), ord = list(o ~ z, family = cumulative()), catg = list(cat ~ z, family = categorical()),
  sigma = list(bf(y ~ z, sigma ~ x)), mv = list(bf(y ~ z) + bf(y2 ~ x) + set_rescor(TRUE)), ar = list(y ~ z + ar(time = t, gr = id)), mo = list(y ~ mo(o) + z),
  me = list(y ~ me(zobs, sdz)), mi = list(bf(y ~ mi(zm)) + bf(zm | mi() ~ 1) + set_rescor(FALSE)), mix = list(bf(y ~ 1), family = mixture(gaussian(), gaussian())),
  nl = list(bf(y ~ a + b * z, a + b ~ 1, nl = TRUE)), hsgp = list(y ~ gp(x, k = 10)), hsgp2 = list(y ~ gp(x, z, k = 6)), hsgpby = list(y ~ gp(x, k = 10, by = g)), hsgpp = list(cnt ~ gp(x, k = 10), family = poisson()),
  reml = list(y ~ z + (1 | g), REML = TRUE), prof = list(y ~ z + (1 | g), control = frmtmb_control(profile = TRUE)), hsgpre = list(y ~ gp(x, k = 10) + (1 | g)), cs = list(y ~ cs(0 + g | h))
)
out <- character()
add <- function(lab, v) { v <- suppressWarnings(as.numeric(unlist(v))); out <<- c(out, sprintf("%s %s", lab, paste(sprintf("%a", v), collapse = ","))) }
for (nm in names(sp)) {
  s <- sp[[nm]]
  f <- tryCatch(suppressWarnings(suppressMessages(do.call(frm, c(s[1], list(data = d), s[-1])))), error = function(e) e)
  if (inherits(f, "error")) { out <- c(out, sprintf("%s ERROR %s", nm, substr(conditionMessage(f), 1, 80))); next }
  add(paste(nm, "opt"), f$opt$par); add(paste(nm, "obj"), f$opt$objective); add(paste(nm, "iter"), f$opt$iterations)
  add(paste(nm, "est"), f$estimates[sort(names(f$estimates))])
  add(paste(nm, "se"), tryCatch(sqrt(diag(as.matrix(vcov(f, full = TRUE)))), error = function(e) NA))
  add(paste(nm, "fit"), tryCatch(suppressWarnings(fitted(f)), error = function(e) NA))
  add(paste(nm, "nodes"), tryCatch(c(nrow(RTMB::GetTape(f$obj, "ADFun")$data.frame())), error = function(e) NA))
  add(paste(nm, "par_names"), match(names(f$obj$env$parameters), c("beta", "betad", "b", "theta", "thetaac", "thetar", "miss")))
}
writeLines(out, file.path(REV_OUT, sprintf("regress-%s-%s.txt", arm, commandArgs(TRUE)[2])))
cat(sprintf("REG %s rows %d errors %d\n", arm, length(out), sum(grepl(" ERROR ", out))))
print(grep(" ERROR ", out, value = TRUE))
