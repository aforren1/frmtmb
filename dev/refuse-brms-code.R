# What brms 2.23.0 writes for the swept pairs frmtmb refuses: the model
# block's density lines. Seed 4040; stancode() needs no compile.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(4040)
n <- 60
d <- data.frame(x = rnorm(n), w = rnorm(n), g = factor(rep(1:6, 10)))
d$y <- rnorm(n); d$y2 <- rnorm(n); d$ym <- d$y; d$ym[1:5] <- NA
d$ce <- sample(c(0L, 0L, 1L), n, TRUE); d$s <- runif(n, .2, .5)
d$wt <- runif(n, .5, 2); d$yo <- sample(1:3, n, TRUE)
d$yp <- exp(d$y)
show <- function(lab, ...) {
  cat("\n=====", lab, "\n")
  code <- tryCatch(suppressWarnings(as.character(brms::stancode(...))),
                   error = function(e) conditionMessage(e))
  ln <- strsplit(code, "\n")[[1]]
  i <- grep("^model \\{", ln)
  if (!length(i)) return(cat(code, "\n"))
  p <- grep("^parameters \\{", ln)
  tp <- grep("^transformed parameters \\{", ln)
  cat(ln[p:(tp - 1L)], sep = "\n")
  cat(ln[i:length(ln)], sep = "\n")
}
mix <- mixture(gaussian(), gaussian())
show("cens mix", bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), family = mix, data = d)
show("trunc mix", bf(y | trunc(lb = -5) ~ 1, mu1 ~ x, mu2 ~ x), family = mix,
     data = d)
show("mi cens", bf(ym | mi() + cens(ce) ~ x), data = d)
show("mi trunc", bf(ym | mi() + trunc(lb = -5) ~ x), data = d)
show("mi se", bf(ym | mi() + se(s) ~ x), data = d)
show("mi(s) se", bf(ym | mi(s) + se(s) ~ x), data = d)
show("skew se", bf(y | se(s) ~ x), family = skew_normal(), data = d)
show("skew se sigma", bf(y | se(s, sigma = TRUE) ~ x), family = skew_normal(),
     data = d)
show("rescor w", bf(y | weights(wt) ~ x) + bf(y2 ~ x) + set_rescor(TRUE),
     data = d)
show("rescor w both", bf(y | weights(wt) ~ x) + bf(y2 | weights(s) ~ x) +
       set_rescor(TRUE), data = d)
show("rescor se", bf(y | se(s) ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = d)
show("rescor mi", bf(ym | mi() ~ x) + bf(y ~ x) + set_rescor(TRUE), data = d)
show("thres mv", bf(yo | thres(2) ~ x, family = cumulative()) + bf(y ~ x) +
       set_rescor(FALSE), data = d)
show("id cor F", bf(y ~ x + (1 + w | gr(g, id = "p", cor = FALSE)),
                    sigma ~ (1 | gr(g, id = "p", cor = FALSE))), data = d)
show("id mix cor", bf(y ~ x + (1 + w | gr(g, id = "p", cor = FALSE)),
                      sigma ~ (1 | gr(g, id = "p"))), data = d)
show("mm pw", bf(y ~ x + (1 | mm(g, g, pw = cbind(wt, wt)))), data = d)
